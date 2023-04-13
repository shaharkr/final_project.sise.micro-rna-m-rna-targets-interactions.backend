import os
import time
def auto_delete():
    path = "C:\\datasets\\temp_files"
    current_time = time.time()
    for file_name in os.listdir(path):
        file_path = os.path.join(path, file_name)
        if os.path.isfile(file_path) and current_time - os.path.getmtime(file_path) > 3600:
            os.remove(file_path)
auto_delete()