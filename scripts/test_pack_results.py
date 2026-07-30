from pathlib import Path


def test_packer_defaults_target_csv_before_using_it() -> None:
    source = (
        Path(__file__).resolve().parent.parent
        / "cloud"
        / "recovery_142"
        / "pack_results.sh"
    ).read_text(encoding="utf-8")
    assignment = 'TARGET_CSV="${TARGET_CSV:-targets_missing_launchable.csv}"'
    assert assignment in source
    assert source.index(assignment) < source.index('[ -e "$PWD/$TARGET_CSV" ]')


def test_packer_activates_alderaan_for_noninteractive_packaging() -> None:
    source = (
        Path(__file__).resolve().parent.parent
        / "cloud"
        / "recovery_142"
        / "pack_results.sh"
    ).read_text(encoding="utf-8")
    assert 'if [ "${CONDA_DEFAULT_ENV:-}" != "alderaan" ]; then' in source
    assert 'conda activate alderaan' in source
