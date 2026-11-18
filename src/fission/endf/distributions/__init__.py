"""Secondary-particle data: how many neutrons come out of a reaction, and with
what energies and directions.

    nu.py          MF1 MT452/455/456  average neutrons per fission, delayed groups,
                   MF6 MT18           multiplicity distribution P(nu, E)
    angular.py     MF4                angular distributions f(mu, E)
    energy.py      MF5                energy distributions f(E -> E')
    mf6.py         MF6                energy-angle distributions (neutron products)
    kinematics.py  two-body kinematics (CM <-> LAB), Kalbach-Mann systematics
    legendre.py    Legendre polynomials and moments
    incident.py    interpolation between incident energies (incl. unit base)

This is a subpackage on purpose: the pointwise-table cache fingerprints the
modules of fission.endf itself, and these readers do not affect the tables.
"""
