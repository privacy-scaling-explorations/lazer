from labrador import poly_t
from fhew import Gadget
from treethings import decompose

class GetMostSignificantBit(Gadget):
    def __init__(self, input: poly_t, output: poly_t, num_bits: int):
        self.input = input
        self.output = output
        self.num_bits = num_bits
        self.bits = []

    def apply(self, builder):
        stat_left = [builder.one]
        for i in range(self.num_bits):
            stat_left += [-builder.constant(2 ** i)]
        # (in, bits[0], ..., bits[n−1])
        witness = [self.input] + self.bits
        builder.PS.fresh_statement(stat_left, witness, builder.zero)

    def compute_witness(self, builder):
        self.bits = [builder.zero] * self.num_bits
        decomposed = decompose(self.input, 2, self.num_bits)

