from pathlib import Path
from typing import Optional
import typer

from vast_post_processing.neighbours import (
    read_release_epochs,
    find_vast_neighbours_by_release_epoch,
)


app = typer.Typer()


@app.command()
def main(
    release_epoch: str,
    vast_data_root: Path = typer.Argument(
        ...,
        help=(
            "Path to VAST data. Must follow the VAST data organization and naming"
            " scheme."
        ),
        exists=True,
        file_okay=False,
        dir_okay=True,
    ),
    release_epochs_csv: Path = typer.Argument(
        ...,
        help=(
            "Path to CSV file containing the release epoch for each VAST image. Each"
            " row must contain at least the observation epoch, field name, SBID, and"
            " release epoch."
        ),
        exists=True,
        file_okay=True,
        dir_okay=False,
    ),
    output_root: Path = typer.Argument(
        ...,
        help="Directory to write output links organized by release epoch and field.",
    ),
    vast_db_repo: Path = typer.Argument(
        ...,
        help="Path to VAST ASKAP Surveys database repository.",
        exists=True,
        file_okay=False,
        dir_okay=True,
    ),
    racs_db_repo: Optional[Path] = typer.Argument(
        None,
        help="Path to RACS ASKAP Surveys database repository.",
        exists=True,
        file_okay=False,
        dir_okay=True,
    ),
    overlap_frac_thresh: float = typer.Option(
        0.05,
        help=(
            "Exclude fields that overlap by less than the given fractional overlap."
            " e.g. the default of 0.05 will exclude fields that overlap by less than 5%"
            " of the area of the central field."
        ),
    ),
    use_corrected: bool = typer.Option(
        True, help="Use the corrected versions of images."
    ),
    neighbours_output: Optional[Path] = typer.Option(
        None,
        help="Write the fields and their neighbours out to a CSV file at the given path.",
        file_okay=True,
        dir_okay=False,
        writable=True,
    ),
    make_links: bool = typer.Option(
        True,
        help=(
            "Make symlinks to images in `output_root`. Default is to make the links."
            " Turn it off with --no-make-links."
        ),
    ),
):
    # get the release epochs
    release_epochs = read_release_epochs(release_epochs_csv)
    # get the neighbours DataFrame and filter for the requested release epoch and
    # overlap area threshold
    vast_neighbours_df = find_vast_neighbours_by_release_epoch(
        release_epoch,
        vast_data_root,
        vast_db_repo,
        release_epochs,
        racs_db_repo=racs_db_repo,
        use_corrected=use_corrected,
    ).query(
        "release_epoch_a == @release_epoch and overlap_frac >= @overlap_frac_thresh"
    )

    if neighbours_output is not None:
        vast_neighbours_df[
            [
                "field_a",
                "sbid_a",
                "obs_epoch_a",
                "release_epoch_a",
                "field_b",
                "sbid_b",
                "obs_epoch_b",
                "release_epoch_b",
                "overlap_frac",
                "delta_t_days",
            ]
        ].to_csv(neighbours_output, index=False)

    # create a directory for each field and create links to the neighbouring images
    if make_links:
        release_output_path = output_root / release_epoch
        release_output_path.mkdir(parents=True, exist_ok=True)
        for _, obs_pair in vast_neighbours_df.iterrows():
            # create directories
            field_inputs_path_a = release_output_path / obs_pair.field_a / "inputs"
            field_inputs_path_a.mkdir(parents=True, exist_ok=True)
            field_inputs_path_b = release_output_path / obs_pair.field_b / "inputs"
            field_inputs_path_b.mkdir(parents=True, exist_ok=True)

            # create a hard link for each field in the pair in both directions, e.g.
            # A/inputs/A.fits, A/inputs/B.fits, B/inputs/A.fits, B/inputs/B.fits (plus weights)
            for output_path in (field_inputs_path_a, field_inputs_path_b):
                target_image_a = output_path / obs_pair.image_path_a.name
                target_weights_a = output_path / obs_pair.weights_path_a.name
                if not target_image_a.exists():
                    obs_pair.image_path_a.link_to(target_image_a)
                if not target_weights_a.exists():
                    obs_pair.weights_path_a.link_to(target_weights_a)

                target_image_b = output_path / obs_pair.image_path_b.name
                target_weights_b = output_path / obs_pair.weights_path_b.name
                if not target_image_b.exists():
                    obs_pair.image_path_b.link_to(target_image_b)
                if not target_weights_b.exists():
                    obs_pair.weights_path_b.link_to(target_weights_b)
