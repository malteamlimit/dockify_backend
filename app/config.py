from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    max_cst_penalty: float = 15
    max_qed_penalty: float = 0.4


settings = Settings()
