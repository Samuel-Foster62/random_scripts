#!/bin/bash
set -e
# From https://zenodo.org/records/7462145 - might contain some useful nuggets for processing
# Processing script for hemozoin datasets collected on the DECTRIS
# SINGLA detector on loan to eBIC in 2022. This version operates on
# .nxs files, which have been produced in NeXus format by Nexgen.

# Check script input
if [ "$#" -ne 1 ]; then
    echo "You must supply the location of the data directory only"
    exit 1
fi

# Set up directories
PROCDIR=$(pwd)
DATADIR=$(realpath "$1")
if [ ! -d "$DATADIR" ]; then
    echo "$DATADIR is not found"
    exit 1
fi
if [[ "$PROCDIR" -ef "$DATADIR" ]]; then
    echo "Please process in a new location, not the data directory"
    exit
fi

cd "$PROCDIR"

# Unit cell restraint
cat > restraint.phil <<+
refinement
{
  parameterisation
  {
    crystal
    {
      unit_cell
      {
        restraints
        {
          tie_to_target
          {
            values=12.086,14.6216,7.9942,90.758,97.093,97.060
            sigmas=0.05,0.05,0.05,0.05,0.05,0.05
          }
        }
      }
    }
  }
}
+


################################
# Generic processing functions #
################################

process_one_singla () {

    DATASET=$1
    D_MIN=$2
    MAX_LATTICES=$3

    echo
    echo "*** PROCESSING IN $(pwd) ***"

    set -x
    dials.import "$DATASET" > /dev/null
    dials.background imported.expt output.plot=background.png\
      n_checkpoints=4 d_max=4 d_min=2 > /dev/null
    dials.find_spots imported.expt d_max=9 d_min="$D_MIN" > /dev/null
    dials.index strong.refl imported.expt detector.fix=distance\
      unit_cell=7.9942,12.086,14.6216,97.060,90.758,97.093\
      max_lattices="$MAX_LATTICES" > /dev/null
    if [ ! -f indexed.expt ]; then
    return 0
    fi
    dials.reindex indexed.expt indexed.refl\
      change_of_basis_op=k,l,h > /dev/null
    dials.refine reindexed.expt reindexed.refl scan_varying=False\
        "$PROCDIR"/restraint.phil  > /dev/null
    dials.show refined.expt | grep "distance" > distance.txt
    dials.show refined.expt | grep "px:" > refined_beam_centre.txt
    dials.refine refined.expt refined.refl detector.fix=distance\
      crystal.unit_cell.force_static=true  > /dev/null
    dials.plot_scan_varying_model refined.expt > /dev/null
    dials.integrate refined.expt refined.refl\
      prediction.d_min="$D_MIN" > /dev/null
    { set +x; } 2>/dev/null
}


integrate () {

    cd "$PROCDIR"

    # Pos2
    # Spots mainly at the last third of the scan. 3 lattices index, but
    # take just the main one
    NAME="20220822_1517_nav14_Pos02"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 1
    cd "$PROCDIR"

    # Pos5
    # Single crystal
    NAME="20220823_1257_nav17_Pos05"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 1
    cd "$PROCDIR"

    # Pos6
    # Single crystal
    NAME="20220823_1259_nav18_Pos06"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 1
    cd "$PROCDIR"

    # Pos8
    # 2 lattices index
    NAME="20220823_1303_nav20_Pos08"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 2
    cd "$PROCDIR"

    # Pos10
    # 3 lattices index
    NAME="20220823_1307_nav22_Pos10"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 3
    cd "$PROCDIR"

    # Pos11
    # Single crystal
    NAME="20220823_1310_nav23_Pos11"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 1
    cd "$PROCDIR"

    # Pos12
    # Single crystal
    NAME="20220823_1311_nav24_Pos12"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 1
    cd "$PROCDIR"

    # Pos13
    # 2 lattices index
    NAME="20220823_1313_nav34_Pos13"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 2
    cd "$PROCDIR"

    # Pos14
    # Single crystal plus minor powder rings esp at 4.1 and 3.6 Å
    NAME="20220823_1315_nav35_Pos14"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 1
    cd "$PROCDIR"

    # Pos16
    # 3 lattices index
    NAME="20220823_1318_nav37_Pos16"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 3
    cd "$PROCDIR"

    # Pos17
    # 2 lattices index
    NAME="20220823_1320_nav38_Pos17"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 2
    cd "$PROCDIR"

    # Pos18
    # 3 lattices index
    NAME="20220823_1322_nav39_Pos18"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 3
    cd "$PROCDIR"

    # Pos19
    # 2 lattices index
    NAME="20220823_1324_nav40_Pos19"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 2
    cd "$PROCDIR"

    # Pos20
    # 4 lattices index
    NAME="20220823_1325_nav41_Pos20"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 4
    cd "$PROCDIR"

    # Pos21
    # 3 lattices plus minor powder rings esp at 4.1 and 3.6 Å
    NAME="20220823_1327_nav42_Pos21"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 3
    cd "$PROCDIR"

    # Pos23
    # Single crystal. First ~145 images blank
    NAME="20220823_1333_nav50_Pos23"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 1
    cd "$PROCDIR"

    # Pos26
    # Single crystal. Strong.
    NAME="20220823_1337_nav53_Pos26"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 1
    cd "$PROCDIR"

    # Pos27
    # 2 lattices index
    NAME="20220823_1339_nav54_Pos27"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 2
    cd "$PROCDIR"

    # Pos28
    # 3 lattices index
    NAME="20220823_1341_nav55_Pos28"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 3
    cd "$PROCDIR"

    # Pos30
    # 2 lattices index
    NAME="20220823_1345_nav57_Pos30"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 2
    cd "$PROCDIR"

    # Pos31
    # 2 lattices index
    NAME="20220823_1346_nav58_Pos31"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 2
    cd "$PROCDIR"

    # Pos32
    # Single lattice indexes.
    NAME="20220823_1348_nav59_Pos32"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 1
    cd "$PROCDIR"

    # Pos33
    # 3 lattices index
    NAME="20220823_1350_nav60_Pos33"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 3
    cd "$PROCDIR"

    # Pos34
    # 7 lattices found!
    NAME="20220823_1352_nav61_Pos34"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 7
    cd "$PROCDIR"

    # Pos35
    # 2 lattices found
    NAME="20220823_1355_nav62_Pos35"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 2
    cd "$PROCDIR"

    # Pos37
    # Single crystal. Fairly weak.
    NAME="20220823_1359_nav64_Pos37"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 1
    cd "$PROCDIR"

    # Pos38
    # 5 lattices index
    NAME="20220823_1401_nav65_Pos38"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 5
    cd "$PROCDIR"

    # Pos39
    # Single crystal
    NAME="20220823_1403_nav66_Pos39"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 1
    cd "$PROCDIR"

    # Pos40
    # 2 lattices index
    NAME="20220823_1404_nav67_Pos40"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 2
    cd "$PROCDIR"

    # Pos41
    # 3 lattices index, but the latter two have few reflections, so
    # here take only the first lattice
    NAME="20220823_1406_nav68_Pos41"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 1
    cd "$PROCDIR"

    # Pos43
    # Single crystal with powder ring esp at 4.1 and 3.6 Å. Weak.
    NAME="20220823_1409_nav70_Pos43"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 1
    cd "$PROCDIR"

    # Pos44
    # Single crystal
    NAME="20220823_1411_nav71_Pos44"
    mkdir -p "$NAME" && cd "$NAME"
    DATASET=$DATADIR/"$NAME".nxs
    process_one_singla "$DATASET" 0.7 1
    cd "$PROCDIR"

}


joint_scale () {

    cd "$PROCDIR"

    set -x

    mkdir -p joint_scale
    cd joint_scale

    dials.combine_experiments\
        "$PROCDIR"/20220822_1517_nav14_Pos02/integrated.expt\
        "$PROCDIR"/20220822_1517_nav14_Pos02/integrated.refl\
        "$PROCDIR"/20220823_1257_nav17_Pos05/integrated.expt\
        "$PROCDIR"/20220823_1257_nav17_Pos05/integrated.refl\
        "$PROCDIR"/20220823_1259_nav18_Pos06/integrated.expt\
        "$PROCDIR"/20220823_1259_nav18_Pos06/integrated.refl\
        "$PROCDIR"/20220823_1303_nav20_Pos08/integrated.expt\
        "$PROCDIR"/20220823_1303_nav20_Pos08/integrated.refl\
        "$PROCDIR"/20220823_1307_nav22_Pos10/integrated.expt\
        "$PROCDIR"/20220823_1307_nav22_Pos10/integrated.refl\
        "$PROCDIR"/20220823_1310_nav23_Pos11/integrated.expt\
        "$PROCDIR"/20220823_1310_nav23_Pos11/integrated.refl\
        "$PROCDIR"/20220823_1311_nav24_Pos12/integrated.expt\
        "$PROCDIR"/20220823_1311_nav24_Pos12/integrated.refl\
        "$PROCDIR"/20220823_1313_nav34_Pos13/integrated.expt\
        "$PROCDIR"/20220823_1313_nav34_Pos13/integrated.refl\
        "$PROCDIR"/20220823_1315_nav35_Pos14/integrated.expt\
        "$PROCDIR"/20220823_1315_nav35_Pos14/integrated.refl\
        "$PROCDIR"/20220823_1318_nav37_Pos16/integrated.expt\
        "$PROCDIR"/20220823_1318_nav37_Pos16/integrated.refl\
        "$PROCDIR"/20220823_1320_nav38_Pos17/integrated.expt\
        "$PROCDIR"/20220823_1320_nav38_Pos17/integrated.refl\
        "$PROCDIR"/20220823_1322_nav39_Pos18/integrated.expt\
        "$PROCDIR"/20220823_1322_nav39_Pos18/integrated.refl\
        "$PROCDIR"/20220823_1324_nav40_Pos19/integrated.expt\
        "$PROCDIR"/20220823_1324_nav40_Pos19/integrated.refl\
        "$PROCDIR"/20220823_1325_nav41_Pos20/integrated.expt\
        "$PROCDIR"/20220823_1325_nav41_Pos20/integrated.refl\
        "$PROCDIR"/20220823_1327_nav42_Pos21/integrated.expt\
        "$PROCDIR"/20220823_1327_nav42_Pos21/integrated.refl\
        "$PROCDIR"/20220823_1333_nav50_Pos23/integrated.expt\
        "$PROCDIR"/20220823_1333_nav50_Pos23/integrated.refl\
        "$PROCDIR"/20220823_1337_nav53_Pos26/integrated.expt\
        "$PROCDIR"/20220823_1337_nav53_Pos26/integrated.refl\
        "$PROCDIR"/20220823_1339_nav54_Pos27/integrated.expt\
        "$PROCDIR"/20220823_1339_nav54_Pos27/integrated.refl\
        "$PROCDIR"/20220823_1341_nav55_Pos28/integrated.expt\
        "$PROCDIR"/20220823_1341_nav55_Pos28/integrated.refl\
        "$PROCDIR"/20220823_1345_nav57_Pos30/integrated.expt\
        "$PROCDIR"/20220823_1345_nav57_Pos30/integrated.refl\
        "$PROCDIR"/20220823_1346_nav58_Pos31/integrated.expt\
        "$PROCDIR"/20220823_1346_nav58_Pos31/integrated.refl\
        "$PROCDIR"/20220823_1348_nav59_Pos32/integrated.expt\
        "$PROCDIR"/20220823_1348_nav59_Pos32/integrated.refl\
        "$PROCDIR"/20220823_1350_nav60_Pos33/integrated.expt\
        "$PROCDIR"/20220823_1350_nav60_Pos33/integrated.refl\
        "$PROCDIR"/20220823_1352_nav61_Pos34/integrated.expt\
        "$PROCDIR"/20220823_1352_nav61_Pos34/integrated.refl\
        "$PROCDIR"/20220823_1355_nav62_Pos35/integrated.expt\
        "$PROCDIR"/20220823_1355_nav62_Pos35/integrated.refl\
        "$PROCDIR"/20220823_1359_nav64_Pos37/integrated.expt\
        "$PROCDIR"/20220823_1359_nav64_Pos37/integrated.refl\
        "$PROCDIR"/20220823_1401_nav65_Pos38/integrated.expt\
        "$PROCDIR"/20220823_1401_nav65_Pos38/integrated.refl\
        "$PROCDIR"/20220823_1403_nav66_Pos39/integrated.expt\
        "$PROCDIR"/20220823_1403_nav66_Pos39/integrated.refl\
        "$PROCDIR"/20220823_1404_nav67_Pos40/integrated.expt\
        "$PROCDIR"/20220823_1404_nav67_Pos40/integrated.refl\
        "$PROCDIR"/20220823_1406_nav68_Pos41/integrated.expt\
        "$PROCDIR"/20220823_1406_nav68_Pos41/integrated.refl\
        "$PROCDIR"/20220823_1409_nav70_Pos43/integrated.expt\
        "$PROCDIR"/20220823_1409_nav70_Pos43/integrated.refl\
        "$PROCDIR"/20220823_1411_nav71_Pos44/integrated.expt\
        "$PROCDIR"/20220823_1411_nav71_Pos44/integrated.refl

    # Workaround "Duplicate batch offsets detected" error
    # https://github.com/xia2/xia2/issues/430
    dials.split_experiments combined.{expt,refl}
    dials.combine_experiments split_*.{expt,refl}

    dials.cosym combined.expt combined.refl

    # Put the cell in the reference setting
    dials.reindex symmetrized.expt symmetrized.refl change_of_basis_op=k,l,h

    # Scale with automated image exclusions, absorption forced off
    dials.scale reindexed.expt reindexed.refl d_min=0.71\
        filtering.method=deltacchalf\
        deltacchalf.mode=image_group\
        deltacchalf.group_size=20\
        deltacchalf.stdcutoff=2.0\
        deltacchalf.max_percent_removed=50\
        physical.absorption_correction=False\
        best_unit_cell=12.086,14.6216,7.9942,90.758,97.093,97.060\
        output.experiments=scaled.expt output.reflections=scaled.refl\
        output.log=scaled.log output.html=scaled.html\
        unmerged_mtz=scaled.mtz merged_mtz=merged.mtz

    # Calculate R_{Friedel} (optional, because dev.dials.r_friedel is
    # not in the path for DIALS release bundles)
    if command -v dev.dials.r_friedel &> /dev/null
    then
        dev.dials.r_friedel hklin=scaled.mtz
    fi

    # Export HKL file for SHELX
    dials.export scaled.{expt,refl} format=shelx shelx.hklout=dials.hkl

    cd "$PROCDIR"
    { set +x; } 2>/dev/null

}


# Image processing
integrate
joint_scale
