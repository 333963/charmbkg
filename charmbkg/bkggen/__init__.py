"""charmbkg.bkggen — generatory i próbki kontrolne tła pod pikiem D (scalenie AA: S + S2 + S4 + S6).

Moduły (bez zmian treści, tylko zebrane):
  S : gen (BkgGenerator/CopulaGen/FlowGen/MassDensity), control, data, study, validate, cli
  S2: flow (FlowBkg/CopulaBkg/MassModel/derive_H/cathode_scores), kin, prep
  S4: feta (FetaTransport, SaladReweighter)
  S6: candrot, kin3
Nowe (AA): api (jeden interfejs Generator.fit(SB)/sample(m)), c2st (jeden test dwóch próbek).
Import leniwy — torch/zuko tylko gdy użyte.
"""
from .api import Generator, available, make_generator  # noqa: F401
from .c2st import c2st  # noqa: F401
