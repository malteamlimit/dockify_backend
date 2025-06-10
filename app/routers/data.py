# import json
#
# from fastapi import APIRouter, Depends
# from sqlmodel import Session, select
#
# from ..db.db import get_session
# from ..models import *
# from ..models_db import *
#
# router = APIRouter()
# SessionDep = Depends(get_session)
#
#
# @router.get("/data", tags=['data'])
# async def get_data(session: Session = Depends(get_session)) -> list[DockingJobPreview]:
#     res = []
#     statement = select(DockingJobRecord)
#     results = session.exec(statement).all()
#
#     for record in results:
#         request = DockingJobRequest.model_construct(**json.loads(record.request))
#         result = DockingResult.model_construct(**json.loads(record.result))
#         complex_list = []
#         for complex_record in result.complex_results:
#             complex_list.append(DockingResultComplex.model_construct(**complex_record))
#         result.complex_results = complex_list
#         pre = DockingJobPreview(
#             job_id=record.job_id,
#             created=request.created,
#             runs=request.runs,
#             best_complex=result.complex_results[result.best_complex_id],
#             weight=result.ligand_properties['weight'],
#             hbond_acc=result.ligand_properties['hbond_acc'],
#             hbond_don=result.ligand_properties['hbond_don'],
#             logp=result.ligand_properties['logp'],
#             qed=result.ligand_properties['qed'],
#         )
#         res.append(pre)
#     return res
#
#
# @router.get("/data/{job_id}", tags=['data'])
# async def get_job_by_id(job_id: str, session: Session = Depends(get_session)) -> DockingJobPublic:
#     statement = select(DockingJobRecord).where(DockingJobRecord.job_id == job_id)
#     result = session.exec(statement).first()
#     if not result:
#         return {"error": "Job not found"}
#
#     request = DockingJobRequest.model_construct(**json.loads(result.request))
#     result_data = DockingResult.model_construct(**json.loads(result.result))
#     request.created = datetime.fromisoformat(request.created)
#     complex_list = []
#     for complex_record in result_data.complex_results:
#         temp = DockingResultComplex.model_construct(**complex_record)
#         temp.violation = deserialize_sets(temp.violation)
#         complex_list.append(temp)
#     result_data.complex_results = complex_list
#     result_data.ligand_properties = DockingResultLigand.model_construct(**result_data.ligand_properties)
#
#     return DockingJobPublic(
#         job_id=result.job_id,
#         request=request,
#         result=result_data
#     )
#
#
# @router.delete("/data/{job_id}", tags=['data'])
# async def delete_job_by_id(job_id: str, session: Session = Depends(get_session)):
#     statement = select(DockingJobRecord).where(DockingJobRecord.job_id == job_id)
#     result = session.exec(statement).first()
#     if not result:
#         return {"error": "Job not found"}
#
#     session.delete(result)
#     session.commit()
#     return {"message": "Job deleted successfully"}
#
#
# async def save_job(job: DockingJob, session: Session):
#     record = DockingJobRecord(
#         job_id=job.job_id,
#         request=job.request.model_dump_json(),
#         result=job.result.model_dump_json()
#     )
#     session.add(record)
#     session.commit()
#
#
# def serialize_sets(obj):
#     if isinstance(obj, set):
#         return list(obj)
#     return obj
#
#
# def deserialize_sets(obj):
#     if isinstance(obj, list):
#         return set(obj)
#     return obj
#
