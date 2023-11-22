import os

def rename_files(folder_path):
    # Get the list of files in the folder
    files = os.listdir(folder_path)

    for file_name in files:
        # Check if the file name contains "_1s"
        if "_1s" in file_name:
            # Construct the new file name by removing "_1s"
            new_name = file_name.replace("_1s", "")

            # Construct the full paths for the old and new names
            old_path = os.path.join(folder_path, file_name)
            new_path = os.path.join(folder_path, new_name)

            # Rename the file
            os.rename(old_path, new_path)

            print(f'Renamed: {file_name} -> {new_name}')

# Specify the path to the folder containing the files
folder_path = "Lab Craters/November 2023 Results II/2023_08_07_50mTorr_h10_1s_860gs_crater16"

# Call the function to rename files in the specified folder
rename_files(folder_path)
