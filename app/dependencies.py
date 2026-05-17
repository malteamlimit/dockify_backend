from . import docking
import threading

thread_local_data = threading.local()

def get_docking_wrapper():
    if not hasattr(thread_local_data, "docking_wrapper"):
        print(30 * "*", "Creating new docking wrapper instance.")
        thread_local_data.docking_wrapper = docking.DockingWrapper()
    return thread_local_data.docking_wrapper
