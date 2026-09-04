import subprocess
import argparse

def wrap_and_center_protein(input_trajectory, output_trajectory, topology, group_number=1):
    """
    Center the protein and apply molecule-level periodic boundary wrapping.

    This performs the same operation as:

        gmx trjconv -pbc mol -center -ur compact

    The protein group is selected automatically for both:
        1. the centering group
        2. the output group

    Arguments
        input_trajectory : Path
            Input trajectory file.

        output_trajectory : Path
            Output trajectory file.

        topology : Path
            GROMACS .tpr file.

        group_number : int, optional
            GROMACS index number corresponding to the Protein group.
            Default is 1.

    Returns
        None
    """

    command = [
        "gmx",
        "trjconv",
        "-f", str(input_trajectory),
        "-s", str(topology),
        "-o", str(output_trajectory),
        "-ur", "compact",
        "-pbc", "mol",
        "-center",
    ]

    # trjconv needs:
    #   1. group to center
    #   2. group to write
    selections = f"{group_number}\n{group_number}\n"

    subprocess.run(command, input=selections, text=True, check=True)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Wraps and centers a protein.")

    parser.add_argument("--raw_traj_path", required=True, help="Path to the input data.") 
    parser.add_argument("--processed_traj_path", required=True, help="Path to the input data.") 
    parser.add_argument("--sequences", required=False, default="F66 M66 V66 L66 A66 Y66 I66", help="Sequence names.")

    # F66 M66 V66 L66 A66 Y66 I66

    args = parser.parse_args()

    raw_traj_path = args.raw_traj_path
    processed_traj_path = args.processed_traj_path
    sequences = args.sequences.split()

    for seq in sequences:
        input_trajectory = (f"{raw_traj_path}/{seq}/{seq}_cut.xtc")
        output_trajectory = (f"{processed_traj_path}/{seq}_wrapped_centered.xtc")
        topology = (f"{raw_traj_path}/{seq}/{seq}.tpr")

        wrap_and_center_protein(input_trajectory, output_trajectory, topology)


