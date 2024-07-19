# models.py
from sqlalchemy import Column, Integer, String, Float, DateTime, ForeignKey, ForeignKeyConstraint, UniqueConstraint, \
    Boolean
from db.connector import Base


class Heatmap(Base):
    __tablename__ = 'heatmap'
    id = Column(Integer, primary_key=True, nullable=False, autoincrement=True)
    temp = Column(Integer)
    u_medida = Column(Float)
    descripcion = Column(Float)
