import os
import datetime

def delete_old_files(folder_path):
    # Get the current time
    current_time = datetime.datetime.now()

    # Iterate over files in the folder
    for file_name in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file_name)

        # Get the file creation time
        creation_time = datetime.datetime.fromtimestamp(os.path.getctime(file_path))

        # Check if the file is older than 1 hour
        if (current_time - creation_time).total_seconds() > 3600:
            os.remove(file_path)
folder_path = "C:\\datasets\\temp_files"
delete_old_files(folder_path)
