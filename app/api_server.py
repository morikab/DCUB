from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import Any, Dict
import sys
import os
import uvicorn

# Ensure modules can be imported
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from modules.main import run_modules

app = FastAPI()

# Allow CORS for local development
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class RunModulesRequest(BaseModel):
    user_input_dict: Dict[str, Any]
    should_run_output_module: bool = True

class RunModulesResponse(BaseModel):
    result: Dict[str, Any]

@app.post("/run-modules", response_model=RunModulesResponse)
def run_modules_endpoint(request: RunModulesRequest):
    try:
        result = run_modules(
            user_input_dict=request.user_input_dict,
            should_run_output_module=request.should_run_output_module
        )
        return RunModulesResponse(result=result)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e)) 


if __name__ == "__main__":
    uvicorn.run(app, host="127.0.0.1", port=8000)
