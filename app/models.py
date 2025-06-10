from enum import Enum

from pydantic import ConfigDict, PrivateAttr
from sqlalchemy import JSON, Column
from sqlmodel import SQLModel, Field, Relationship



class JobStatus(str, Enum):
    DRAFT = "draft"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


class DockingJobBase(SQLModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    job_id: str = Field(primary_key=True)
    name: str
    created: str
    constraints: list[list[int | str | list[float] | float]] | None = Field(default=None, sa_column=Column(JSON))
    job_status: JobStatus = Field(default=JobStatus.DRAFT)
    runs: int = Field(default=0, ge=0)
    smiles: str
    sdf: str | None = None

    error: str | None = None
    progress: int = Field(default=0, ge=0, le=100)
    progress_info: str | None = None
    weight: float | None = None
    hbond_acc: int | None = None
    hbond_don: int | None = None
    logp: float | None = None
    qed: float | None = None
    is_sub: bool = Field(default=True)

    best_complex_nr: int | None = None


class DockingJob(DockingJobBase, table=True):
    complexes: list["ComplexResult"] = Relationship(back_populates="job", cascade_delete=True)

    _listeners: list[callable] = PrivateAttr(default_factory=list)

    def get_listeners(self):
        if not hasattr(self, '__pydantic_private__') or self.__pydantic_private__ is None:
            self.__pydantic_private__ = {}
        if '_listeners' not in self.__pydantic_private__:
            self.__pydantic_private__['_listeners'] = []
        return self.__pydantic_private__['_listeners']

    def update(self, **kwargs):
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)

        for listener in self.get_listeners():
            try:
                listener()
            except Exception as e:
                print(f"[update] Listener failed: {e}")


class ComplexResult(SQLModel, table=True):
    violation: str | None

    job_id: str | None = Field(default=None, primary_key=True, foreign_key="dockingjob.job_id", ondelete="CASCADE")
    id: int | None = Field(default=None, primary_key=True)

    total_score: float
    atom_pair_cst: float
    atom_attraction: float
    electrostatic: float
    atom_repulsion: float
    solvation: float
    hbond: float
    delta_g: float
    pairwise_energy: float
    rmsd: float | None = None

    job: DockingJob = Relationship(back_populates="complexes")


class DockingJobWComp(DockingJobBase):
    model_config = ConfigDict(from_attributes=True)

    complexes: list[ComplexResult] = []