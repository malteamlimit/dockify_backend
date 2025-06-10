from datetime import datetime
from enum import Enum

from pydantic import ConfigDict
from sqlmodel import SQLModel, Field, Relationship


# # -------------- DOCKING JOB RESULTS --------------
# class DockingResultLigandDB(SQLModel, table=True):
#     id: int | None = Field(default=None, primary_key=True)
#     weight: float
#     hbond_acc: int
#     hbond_don: int
#     logp: float
#     qed: float
#     preview_path: str | None = None
#
#     job_id: str = Field(default=None, foreign_key="dockingjobdb.job_id")
#     job: "DockingJobDB" = Relationship(back_populates="ligand")
#
#
# class DockingResultComplexDB(SQLModel, table=True):
#     id: int | None = Field(default=None, primary_key=True)
#     total_score: float
#     atom_pair_cst: float
#     atom_attraction: float
#     electrostatic: float
#     atom_repulsion: float
#     solvation: float
#     hbond: float
#     delta_g: float
#     pairwise_energy: float
#     rmsd: float
#     violation: str | None = None  # Als CSV oder JSON
#
#     pose_path: str | None = None
#
#     job_id: str = Field(default=None, foreign_key="dockingjobdb.job_id")
#     job: "DockingJobDB" = Relationship(back_populates="complexes")
#
#
# # -------------- DOCKING JOB REQUEST --------------
# class DockingJobRequestDB(SQLModel, table=True):
#     job_id: str = Field(primary_key=True)
#     created: datetime.datetime
#     smiles: str
#     constraints: str | None = None
#     runs: int
#
#
# # -------------- DOCKING JOB --------------
# class DockingJobDB(SQLModel, table=True):
#     job_id: str = Field(primary_key=True)
#     best_complex_id: int
#
#     request_id: int | None = Field(default=None, foreign_key="dockingjobrequestdb.job_id")
#     request: DockingJobRequestDB | None = Relationship()
#
#     ligand: DockingResultLigandDB = Relationship(back_populates="job")
#     complexes: list[DockingResultComplexDB] = Relationship(back_populates="job")


# class DockingJobRecord(SQLModel, table=True):
#     job_id: str = Field(primary_key=True)
#     request: str
#     result: str

