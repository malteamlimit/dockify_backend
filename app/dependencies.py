from . import docking
import threading

thread_local_data = threading.local()

def get_docking_wrapper():
    from app.main import main_event_loop

    if not hasattr(thread_local_data, "docking_wrapper"):
        print(30 * "*", "Creating new docking wrapper instance.")
        thread_local_data.docking_wrapper = docking.DockingWrapper(loop=main_event_loop)
    return thread_local_data.docking_wrapper
