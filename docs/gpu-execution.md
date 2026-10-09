# NVIDIA GPU execution

BFQ's GPU path is host driver → Docker → BFQ → Snakemake → Apptainer
`--nv` → the workflow image. Enable it explicitly with `BFQ_GPU=1` in the
Docker launcher. This adds `--singularity-args=--nv` to every BFQ Snakemake
invocation, including analysis resume. Unset or `BFQ_GPU=0` keeps CPU-only
operation; it does not make a GPU-dependent scientific workflow run on a CPU.
Accepted boolean spellings also include true/false, yes/no and on/off.

The BFQ base retains its Miniconda environment. It declares NVIDIA driver
capabilities `compute,utility`, registers conventional NVIDIA library mount
directories, and directs Apptainer's `binary path` to system `ldconfig` before
Conda's tools. The host's NVIDIA Container Toolkit injects the matching driver
libraries at **container launch**, including `libcuda.so.1`. They cannot be
validated during an ordinary image build. No host driver, CUDA stub library or
second CUDA toolkit is installed in BFQ: the scientific image supplies its
runtime. See [Apptainer GPU support](https://apptainer.org/docs/user/main/gpu.html)
and [NVIDIA Docker capabilities](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/docker-specialized.html).

## Host and operational launcher

Install a supported NVIDIA driver and NVIDIA Container Toolkit on the Docker
host using the [NVIDIA installation guide](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html).
Configuring the Docker runtime is a host administration step, for example
`sudo nvidia-ctk runtime configure --runtime=docker`, followed by the planned
Docker restart. Record the actual host evidence before accepting the deployment:

```bash
nvidia-smi --query-gpu=name,uuid,driver_version --format=csv
nvidia-ctk --version
nvidia-container-cli --version
docker version
docker info --format '{{json .Runtimes}}'
```

The workflow-owned `docker.config` selects the RAPIDS image. The integration
failure that motivated #148 used `gcfntnu/rapids-scanpy:0.17.0`, containing CUDA
12.9.1 / CUDA runtime 12.9.79. These are observed image versions, **not a completed
host compatibility test**. CUDA 12.x minor-version compatibility has driver and
feature restrictions, including PTX compilation; a reported driver version or
`nvidia-smi` alone cannot establish that this workload works. Consult
[NVIDIA compatibility guidance](https://docs.nvidia.com/deploy/cuda-compatibility/minor-version-compatibility.html)
and run the computation below with the actual image and GPU.

Add these arguments to the actual test and production `docker run` commands,
before the image name, while retaining their existing operational mounts:

```bash
--gpus all -e NVIDIA_DRIVER_CAPABILITIES=compute,utility -e BFQ_GPU=1
```

The test image is `gcfntnu/bfq:dev-test`. Recreate the container with the updated
launcher to acquire driver libraries; setting environment variables in an
already running container does not add Docker GPU passthrough. CPU-only
launchers omit `--gpus` and set `BFQ_GPU=0` or leave it unset.

The external launchers found during development are:

| Role | Mounted script | Required change |
| --- | --- | --- |
| Reduced integration test candidate | `/mnt/gcf-work/bfq-test/scripts/reduced-data-it-start-bfq.sh` | Add the three GPU arguments to the active `docker run` line; use `dev-test`. |
| Production candidate | `/mnt/gcf-work/bfq-prod/scripts/bfq-it-start.sh` | Add the same arguments when promoting a tested production image. |

Confirm these are the launchers actually used on the target host and record
their final paths, checksums/diffs and invocation in the PR. They are maintained
outside this repository. The reduced test script also clears prior test state
and outputs, so its clean-start path is distinct from resuming an existing run.
Do not use it to prepare a retained-workdir resume.

## Build the two integration images

From the issue checkout, use the existing build-and-push script:

```bash
./build-tag-push.sh base base-test
./build-tag-push.sh test dev-test -b gcfntnu/bfq:base-test -w bfq-dev
```

These commands publish the named test tags. The `-b` option changes only the
test image's base; production builds retain their existing base policy. If the
current account already has Docker access, set `BFQ_DOCKER=docker` for both
commands. Record both registry digests, BFQ source commit and the installed
tools/workflow revisions with the integration results.

After integration approval, publish a separately named production base and
update the base reference in `dockerfile-prod` as part of production promotion.
Its current `base-260925` reference does not acquire these changes when
`base-test` is published. Do not promote either reusable test tag as a production
release.

## Smoke through the real execution path

Start the disposable test container with its normal mounts and the three GPU
arguments above. Before starting BFQ, run inside it:

```bash
bfq-gpu-smoke
```

The command resolves `rapids-scanpy` from the installed workflow `docker.config`
(including the existing `GCF_WORKFLOWS_DOCKER_CONFIG` override). It retains a
unique scratch directory with diagnostics, commands, Snakemake logs and results.
It checks Docker-visible devices, driver-library loading and Apptainer library
discovery, then executes a tiny GPU computation through the **same command
builder used by BFQ analysis and resume**. Both paths must initialize CUDA,
allocate GPU memory, compute, synchronize and validate the result. It uses no
flowcell inputs, state, scientific outputs or email. Run it on the intended GPU
host; the normal local development suite mocks these external operations.

A minimal independent launcher for that smoke check is:

```bash
smoke_root=$(mktemp -d)
docker run --rm --privileged --gpus all \
  -e NVIDIA_DRIVER_CAPABILITIES=compute,utility -e BFQ_GPU=1 \
  -e SINGULARITY_BINDPATH=/bfq-tmp -e APPTAINER_BINDPATH=/bfq-tmp \
  -v "$smoke_root:/bfq-tmp" \
  gcfntnu/bfq:dev-test bfq-gpu-smoke
```

This uses an isolated cache and may download the large RAPIDS image. The final
acceptance run must also pass through the actual operational launcher, with its
mounts and environment. Retain the smoke directory for review.

| Failure | Inspect |
| --- | --- |
| No GPU devices / Docker rejects `--gpus` | Host driver, toolkit/runtime configuration, actual Docker launch arguments. |
| `libcuda.so.1` cannot load or is absent from `ldconfig -p` | Docker driver injection and `compute` capability; system loader paths. |
| Apptainer cannot discover driver libraries | `singularity config global --get 'binary path'`, `/sbin/ldconfig -p`, and rebuilt BFQ base. |
| No nested GPU / missing `--nv` | `BFQ_GPU=1` in the BFQ process and the emitted Snakemake/Apptainer commands. |
| Insufficient driver / unsupported PTX / CUDA initialization failure | Actual runtime, host driver and GPU compatibility; retain the full nested-container error. |
| Mount, permission or image-pull failure | Existing Apptainer privileges, bind paths, scratch capacity and registry access. |

## Integration acceptance and issue closure

Keep #148 open until the target-host evidence is complete. Local tests and
published images are preparation for the maintainer's integration tests.

1. Record host GPU/driver/toolkit/Docker versions, both BFQ image digests, the
   workflow commit and actual RAPIDS image identity. Verify the active external
   test launcher changes and retain both successful smoke results.
2. Run the established reduced-data integration suite with complete current
   single-cell defaults. Do not restore the obsolete notebook-era
   `setup_test/singlecell.config`. Aggregation normalization remains `none` for
   the heavily downsampled data.
3. Confirm `preprocess_native_representation` passes the previously failing
   `rmm.reinitialize` step, and inspect subsequent scientific outputs/reports.
4. Verify retained-workdir analysis resume with the corrected GPU launcher:
   `fm rerun RUN_ID --from analysis --resume --dry-run`, then queue without
   `--dry-run`. Resume reuses the copied workflow and configuration. To test
   newly installed workflow code instead, use an ordinary analysis rerun.
5. Run an existing CPU-only workflow without GPU passthrough / with `BFQ_GPU=0`.
6. Record and verify the production launcher change during deliberate promotion.
   Merging this PR does not deploy images or complete the external acceptance.
