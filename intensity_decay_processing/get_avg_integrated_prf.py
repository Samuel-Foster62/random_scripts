#!/bin/env dials.python

from __future__ import annotations

import os
import sys

import iotbx.phil
from cctbx import uctbx
from dxtbx.model import ExperimentType
from dxtbx.model.experiment_list import ExperimentListFactory
from scitbx.math import five_number_summary

import dials.util
from dials.array_family import flex
from dials.util import Sorry, tabulate

help_message = """

Examples::

  dials.show models.expt

  dials.show image_*.cbf

  dials.show observations.refl
"""

phil_scope = iotbx.phil.parse(
    """\
show_scan_varying = False
  .type = bool
  .help = "Whether or not to show the crystal at each scan point."
show_shared_models = False
  .type = bool
  .help = "Show which models are linked to which experiments"
show_all_reflection_data = False
  .type = bool
  .help = "Whether or not to print individual reflections"
show_intensities = False
  .type = bool
show_centroids = False
  .type = bool
show_profile_fit = False
  .type = bool
show_flags = False
  .type = bool
  .help = "Show a summary table of reflection flags"
show_identifiers = False
  .type = bool
  .help = "Show experiment identifiers map if set"
image_statistics{
  show_corrected = False
    .type = bool
    .help = "Show statistics on the distribution of values in each corrected image"
  show_raw = False
    .type = bool
    .help = "Show statistics on the distribution of values in each raw image"
}
max_reflections = None
  .type = int
  .help = "Limit the number of reflections in the output."
""",
    process_includes=True,
)

def run(args=None):
    import dials.util.log

    dials.util.log.print_banner()

    from dials.util.options import (
        ArgumentParser,
        reflections_and_experiments_from_files,
    )

    usage = "dials.show [options] models.expt | image_*.cbf"

    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        read_experiments_from_images=True,
        read_reflections=True,
        check_format=False,
        epilog=help_message,
    )

    params, options = parser.parse_args(args=args, show_diff_phil=True)
    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    if len(experiments) == 0 and len(reflections) == 0:
        parser.print_help()
        exit()

#    if len(experiments):
#        if not all(e.detector for e in experiments):
#            sys.exit("Error: experiment has no detector")
#        if not all(e.beam for e in experiments):
#            sys.exit("Error: experiment has no beam")
#        print(show_experiments(experiments, show_scan_varying=params.show_scan_varying))
#
#        if params.image_statistics.show_raw:
#            show_image_statistics(experiments, "raw")
#
#        if params.image_statistics.show_corrected:
#            show_image_statistics(experiments, "corrected")
#
#        if params.show_shared_models:
#            print()
#            print(model_connectivity(experiments))

    if len(reflections):
        print(
            show_reflections(
                reflections,
                show_intensities=params.show_intensities,
                show_profile_fit=params.show_profile_fit,
                show_centroids=params.show_centroids,
                show_all_reflection_data=params.show_all_reflection_data,
                show_flags=params.show_flags,
                max_reflections=params.max_reflections,
                show_identifiers=params.show_identifiers,
            )
        )

def show_reflections(
    reflections,
    show_intensities=False,
    show_profile_fit=False,
    show_centroids=False,
    show_all_reflection_data=False,
    show_flags=False,
    max_reflections=None,
    show_identifiers=False,
):
    text = []

    for rlist in reflections:
        if "intensity.prf.value" in rlist:
            mean_intensity = flex.mean(rlist["intensity.prf.value"])
            text.append(f"{mean_intensity:.1f}")
    return "\n".join(text)

if __name__ == "__main__":
    run()