from abc import ABC, abstractmethod
from labrador import poly_t
from lazer import polyring_t

class Gadget(ABC):
    @abstractmethod
    def apply(self, builder: "PrincipalRelationBuilder") -> None:
        """
        Given a PrincipalRelationBuilder (with its PS inside),
        add this constraint to it.
        """
        ...

    @abstractmethod
    def compute_witness(self, builder: "PrincipalRelationBuilder"):
        """
        Compute witness values for the gadget
        that satisfies gadget constraints
        """
        ...

class PrincipalRelationBuilder:
    def __init__(self,
        deg_list: list[int],
        num_pols: list[int],
        norm_list: list[int],
        prime_sz: str,
        modulus: int,
        num_constraints: int
    ):
        self.ring = polyring_t(deg_list[0], modulus)   # assume uniform deg
        self.one  = poly_t.from_int(self.ring, 1)
        self.zero = poly_t.zero(self.ring)

        # create the PS object once
        self.PS = proof_statement(
            deg_list,
            num_pols,
            norm_list,
            num_constraints,
            prime_sz
        )

    def constant(self, constant: int):
        return poly_t.from_int(self.ring, constant)

    def add(self, constraint: Gadget):
        constraint.apply(self)

    def finalize(self):
        stmt = self.PS.output_statement()
        assert self.PS.smpl_verify()
        return self.PS.pack_prove()

