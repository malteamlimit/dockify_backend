import asyncio
from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles
from starlette.middleware.cors import CORSMiddleware
from starlette.responses import HTMLResponse

from . import docking
from .db import db
from .routers import util, jobs, data
# from .dependencies import thread_local_data


# TODO:
# - Add qed validation (constraint: 0,4)
# - Add core validation (is substructure of the core?)

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
# app.include_router(data.router)


# @app.get("/")
# def read_root():
#     return {"API": "running"}


@app.get("/")
async def get():
    return HTMLResponse("""
    
<!DOCTYPE html>
<html>
    <head>
        <title>Chat</title>
    </head>
    <body>
        <h1>WebSocket Test</h1>
        <form action="" onsubmit="sendJob(event)">            
            <button>Send</button>
        </form>
        <ul id='messages'>
        </ul>
        <script>
            function sendJob(event) {
                // send http job to start docking
                event.preventDefault();
                fetch("/jobs/create", {
                    method: "POST",
                    headers: {
                        "Content-Type": "application/json"
                    },
                    body: JSON.stringify({
                        "smiles": "CN(C(=O)CN3CC2(CCN(C(=O)c1cccnc1)CC2)C3)c5ccc4COCc4c5",
                        "constraints": [[364,  'HG', [-6.7520, -0.1555, 13.0855], 1.80, 0.125], [ 65, 'OD2', [-7.1638,  5.8368, 16.5862], 3.23, 0.25], [ 65, 'OD2', [-7.5181,  3.1143, 15.5623], 3.25, 0.25], [ 89,  'CB', [-6.0966,  5.3594, 15.7673], 3.70, 0.25], [ 86,  'CD', [-7.1638,  5.8368, 16.5862], 5.11, 0.50]],
                        "runs": 15
                    })
                })
                .then(response => response.json())
                .then(data => {
                    if (data.job_id === undefined) {
                        return;
                    }
                    console.log(data);
                    // Start WebSocket connection
                    connectWebSocket(data.job_id);
                });
            }
            function connectWebSocket(jobId) {
                const socket = new WebSocket(`ws://localhost:8000/jobs/status/${jobId}`);
                socket.onopen = function(event) {
                    console.log("WebSocket connection opened");
                };
                socket.onmessage = function(event) {
                    const message = JSON.parse(event.data);
                    console.log("Message from server: ", message);
                    const messagesList = document.getElementById('messages');
                    const listItem = document.createElement('li');
                    listItem.textContent = JSON.stringify(message);
                    messagesList.appendChild(listItem);
                };
                socket.onclose = function(event) {
                    console.log("WebSocket connection closed");
                };
                socket.onerror = function(error) {
                    console.error("WebSocket error: ", error);
                };
            }
        </script>
    </body>
</html> 

    """)
