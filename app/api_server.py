from fastapi import FastAPI, File, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field
from typing import Any, Dict
import shutil
import sys
import os
import uvicorn

# Ensure modules can be imported
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from modules.main import artifacts_directory, run_modules
from modules.models import UserInput

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
    user_input: UserInput = Field(..., alias="user_input_dict")
    should_run_output_module: bool = True


class RunModulesResponse(BaseModel):
    result: Dict[str, Any]


@app.post("/upload-file")
async def upload_file(file: UploadFile = File(...)):
    """Dev-only: accept a file upload and return its server-side path for use in /run-modules."""
    upload_dir = artifacts_directory / "uploads"
    upload_dir.mkdir(parents=True, exist_ok=True)
    dest = upload_dir / file.filename
    with dest.open("wb") as f:
        shutil.copyfileobj(file.file, f)
    return {"file_path": str(dest)}


@app.post("/run-modules", response_model=RunModulesResponse)
def run_modules_endpoint(request: RunModulesRequest):
    try:
        result = run_modules(
            user_input=request.user_input,
            should_run_output_module=request.should_run_output_module
        )
        if "error_message" in result:
            raise HTTPException(status_code=500, detail=result["error_message"])
        return RunModulesResponse(result=result)
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


if __name__ == "__main__":
    uvicorn.run(app, host="127.0.0.1", port=8000)
