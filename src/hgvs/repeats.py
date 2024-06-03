# -*- coding: utf-8 -*-
""" A class to manage conversion of SequenceVariants to repeat representation"""

import re
from typing import Optional

from hgvs.dataproviders.interface import Interface
from hgvs.pretty.datacompiler import DataCompiler
from hgvs.pretty.models import PrettyConfig, VariantCoords
from hgvs.sequencevariant import SequenceVariant


def count_pattern_occurence(pattern, string):
    """Counts how often a pattern can be found in the string."""
    matches = re.findall(pattern, string)
    return len(matches)


class RepeatAnalyser:

    def __init__(self, hdp: Interface, var_g: SequenceVariant) -> None:

        config = PrettyConfig(hdp, None, None)
        dc = DataCompiler(config)

        fs = dc.get_shuffled_variant(var_g, 0)

        self.is_repeat = False
        self.ref_count = 0
        self.alt_count = 0
        self.ref_str = fs.ref
        self.alt_str = fs.alt

        self.repeat_unit = self._get_repeat_unit(fs)
        if self.repeat_unit is None:
            return

        self.is_repeat = True

        self.ref_count = count_pattern_occurence(self.repeat_unit, fs.ref)
        self.alt_count = count_pattern_occurence(self.repeat_unit, fs.alt)
        self.ref_str = f"{self.repeat_unit}[{self.ref_count}]"
        self.alt_str = f"{self.repeat_unit}[{self.alt_count}]"

    def __repr__(self):
        return f"{self.ref_str}>{self.alt_str}"

    def _get_repeat_unit(self, fs: VariantCoords) -> Optional[str]:
        """Takes fully justified coordiantes and tries to detect a repeat in them."""
        # analyze for repeat:
        if len(fs.ref) == len(fs.alt):
            # seems we cant shuffle. is an SVN or delins
            return None

        if len(fs.alt) > 0:
            if fs.alt in fs.ref:
                if fs.ref.startswith(fs.alt):
                    d = fs.ref[len(fs.alt) :]
                    if count_pattern_occurence(d, fs.ref) > 1:
                        return d
            elif fs.ref in fs.alt:
                if fs.alt.startswith(fs.ref):
                    d = fs.alt[len(fs.ref) :]
                    if count_pattern_occurence(d, fs.ref) > 1:
                        return d

        return None
