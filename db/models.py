# models.py
from sqlalchemy import Column, Integer, String, Float, DateTime, ForeignKey, ForeignKeyConstraint, UniqueConstraint, \
    Boolean
from db.connector import Base


class Heatmap(Base):
    __tablename__ = 'heatmap'
    id = Column(Integer, primary_key=True, nullable=False, autoincrement=True)
    temp = Column(Float)
    mltss = Column(Integer)
    sludge_prod = Column(Float)

    def to_dict(self):
        return {
            'id': self.id,
            'temp': self.temp,
            'mltss': self.mltss,
            'sludge_prod': self.sludge_prod,
        }