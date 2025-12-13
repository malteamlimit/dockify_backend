import asyncio
from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles
from starlette.middleware.cors import CORSMiddleware

from . import docking
from .db import db
from .routers import util, jobs, database
# from .dependencies import thread_local_data


main_event_loop = None

@asynccontextmanager
async def lifespan(app: FastAPI):
    global main_event_loop
    main_event_loop = asyncio.get_running_loop()
    print("Event Loop initialized:", main_event_loop)
    # thread_local_data.docking_wrapper = docking.DockingWrapper()
    db.init_db()
    print(30 * "*", "Docking wrapper & DB initialized.")
    yield


app = FastAPI(lifespan=lifespan)



origins = [
    "http://localhost:3000",
    "http://127.0.0.1:3000",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.mount("/static", StaticFiles(directory="app/static/"), name="static")

app.include_router(util.router)
app.include_router(jobs.router)
app.include_router(database.router)


@app.get("/")
def read_root():
    return {"API": "running"}
