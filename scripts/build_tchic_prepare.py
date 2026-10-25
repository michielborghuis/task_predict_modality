"""Build the scoped Viash 0.9.7 preparation and an isolated Apptainer runtime.

The R image recipe is derived from process_dataset's declared Docker setup;
no packages are installed outside the new container. Existing SIFs are checked
by content hash. Logs and runtime identities remain under the selected target.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

import yaml

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TARGET = "target_tchic_verified_20261023"
PIN = "106edbdeae9775d71bdfb4e2dcca67e187dc843b"
PYTHON_SHA = "c7a010cf9a282565c9401e74de3438ab33587d770bc52112b4981b6091fdd6fb"


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_dependency():
    dep = ROOT / ".dependencies/datasets-106edbde"
    actual = subprocess.check_output(["git", "--no-optional-locks", "-C", str(dep), "rev-parse", "HEAD"], text=True).strip()
    if actual != PIN:
        raise RuntimeError(f"Expected datasets {PIN}, found {actual}")
    subprocess.run(["git", "--no-optional-locks", "-C", str(dep), "diff", "--quiet", "HEAD", "--",
                    "src/normalization/log_cp", "src/api", "_viash.yaml"], check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", default=DEFAULT_TARGET)
    cache = Path("/home/michiel/.apptainer/cache/nextflow_tchic")
    parser.add_argument("--base-r", type=Path, default=cache / "openproblems-base_r-1.img")
    parser.add_argument("--base-python", type=Path, default=cache / "openproblems-base_python-1.img")
    parser.add_argument("--runtime-only", action="store_true")
    args = parser.parse_args()
    target = (ROOT / args.target).resolve()
    runtime = target / "runtime"
    runtime.mkdir(parents=True, exist_ok=True)
    verify_dependency()
    if "0.9.7" not in subprocess.check_output(["viash", "--version"], text=True):
        raise RuntimeError("This build requires Viash 0.9.7")
    config_path = ROOT / "src/data_processors/process_dataset/config.vsh.yaml"
    config = yaml.safe_load(config_path.read_text())
    expected_r = config["info"]["tchic_runtime"]["base_sif_sha256"]
    for path, expected in [(args.base_r, expected_r), (args.base_python, PYTHON_SHA)]:
        if sha256(path) != expected:
            raise RuntimeError(f"Base image checksum mismatch: {path}")

    def run(command, name):
        print("RUN", command, flush=True)
        prior = (runtime / "commands.jsonl").read_text().splitlines() if (runtime / "commands.jsonl").exists() else []
        log_path = runtime / f"{len(prior) + 1:02d}_{name}"
        with log_path.open("w") as log:
            result = subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
        with (runtime / "commands.jsonl").open("a") as record:
            record.write(json.dumps({"command": command, "exit_code": result.returncode, "log": str(log_path)}) + "\n")
        if result.returncode:
            print(log_path.read_text()[-12000:], flush=True)
            raise RuntimeError(f"Command failed ({result.returncode}); see {log_path}")

    image = runtime / "process_dataset.sif"
    lock_path = runtime / "runtime.json"
    if image.exists():
        lock = json.loads(lock_path.read_text())
        if lock["processor_sha256"] != sha256(image) or lock["processor_config_sha256"] != sha256(config_path):
            raise RuntimeError("Runtime/config changed: use a new target; existing images are not overwritten")
    else:
        docker = next(engine for engine in config["engines"] if engine["type"] == "docker")
        setup = "\n".join(item["run"] for item in docker["setup"])
        definition = runtime / "process_dataset.def"
        definition.write_text(
            f"Bootstrap: localimage\nFrom: {args.base_r.resolve()}\n\n%post\n    set -eu\n" +
            "\n".join("    " + line for line in setup.splitlines()) +
            "\n\n%environment\n    export RETICULATE_PYTHON=/opt/tchic-python/bin/python\n"
            "    export RETICULATE_USE_MANAGED_VENV=no\n")
        run(["apptainer", "build", "--fakeroot", str(image), str(definition)], "image_build.log")
        lock = {"datasets_commit": PIN, "python_image": str(args.base_python.resolve()),
                "python_sha256": PYTHON_SHA, "r_base_image": str(args.base_r.resolve()),
                "r_base_sha256": expected_r, "processor_image": str(image),
                "processor_sha256": sha256(image), "processor_config_sha256": sha256(config_path)}
        lock_path.write_text(json.dumps(lock, indent=2) + "\n")
    run(["apptainer", "exec", "--cleanenv", "--no-home", str(image), "/opt/tchic-python/bin/python", "-m", "pip", "freeze"], "python_packages.txt")
    run(["apptainer", "exec", "--cleanenv", "--no-home", str(image), "Rscript", "-e",
         'cat(R.version.string,"\\n"); for(p in c("anndata","reticulate","Matrix")) cat(p,as.character(packageVersion(p)),"\\n"); print(reticulate::py_config())'], "r_python_versions.txt")
    if args.runtime_only:
        return
    repository = ROOT / ".dependencies/datasets-106edbde"
    dependency_target = repository / "target_viash097"
    dependency_config = dependency_target / "nextflow/normalization/log_cp/.config.vsh.yaml"
    if dependency_config.is_file():
        dependency_info = yaml.safe_load(dependency_config.read_text()).get("build_info", {})
        if dependency_info.get("git_commit") != PIN:
            raise RuntimeError(f"Existing dependency target is not built from {PIN}: {dependency_target}")
        for artifact in [
            dependency_target / "executable/normalization/log_cp/log_cp",
            dependency_target / "nextflow/normalization/log_cp/main.nf",
        ]:
            if not artifact.is_file():
                raise RuntimeError(f"Incomplete pinned dependency target: {artifact}")
        print(f"Reusing pinned dependency target {dependency_target}", flush=True)
    elif dependency_target.exists():
        raise RuntimeError(f"Incomplete dependency target exists; refusing to overwrite it: {dependency_target}")
    else:
        run(["bash", "scripts/bootstrap_tchic_datasets.sh", str(dependency_target)], "build_dependency.log")
    # Viash 0.9.7 resolves local repositories from the checkout root, not the built target.
    mod = '.repositories[.name == "tchic_datasets_106edbde"].path := "/' + str(repository.relative_to(ROOT)) + '"'
    run(["viash", "ns", "build", "--query", "^data_processors/(prepare_tchic_input|process_dataset)$",
         "--target", str(target), "--runner", "nextflow", "--engine", "docker", "--config_mod", mod], "build_nextflow.log")
    run(["viash", "ns", "build", "--query", "^workflows/prepare_tchic$",
         "--target", str(target), "--runner", "nextflow", "--engine", "native", "--config_mod", mod], "build_workflow.log")
    # Native wrappers execute INSIDE the selected images, never in host R/Python.
    run(["viash", "ns", "build", "--query", "^data_processors/(prepare_tchic_input|process_dataset)$",
         "--target", str(target), "--runner", "executable", "--engine", "native",
         "--config_mod", '.engines += [{"type":"native"}]'], "build_executable.log")
    sources = ["src/data_processors/prepare_tchic_input/script.py", "src/data_processors/process_dataset/script.R",
               "src/data_processors/process_dataset/config.vsh.yaml", "src/workflows/prepare_tchic/main.nf",
               "src/workflows/prepare_tchic/config.vsh.yaml", "src/workflows/prepare_tchic/runtime.config",
               "scripts/bootstrap_tchic_datasets.sh"]
    lock["source_sha256"] = {path: sha256(ROOT / path) for path in sources}
    artifacts = ["nextflow/workflows/prepare_tchic/main.nf", "nextflow/data_processors/prepare_tchic_input/main.nf",
                 "nextflow/data_processors/process_dataset/main.nf",
                 "../.dependencies/datasets-106edbde/target_viash097/executable/normalization/log_cp/log_cp",
                 "../.dependencies/datasets-106edbde/target_viash097/nextflow/normalization/log_cp/main.nf",
                 "dependencies/local/nextflow/normalization/log_cp/main.nf",
                 "executable/data_processors/prepare_tchic_input/prepare_tchic_input",
                 "executable/data_processors/process_dataset/process_dataset"]
    lock["artifact_sha256"] = {path: sha256(target / path) for path in artifacts}
    bound = yaml.safe_load((target / "dependencies/local/nextflow/normalization/log_cp/.config.vsh.yaml").read_text())
    assert bound["build_info"]["git_commit"] == PIN
    lock_path.write_text(json.dumps(lock, indent=2) + "\n")
    print(f"Built {target}; identities: {lock_path}", flush=True)


if __name__ == "__main__":
    main()
