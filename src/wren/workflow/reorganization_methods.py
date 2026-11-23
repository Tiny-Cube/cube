"""Available post-processing algorithms; importing controls never loads a QC engine."""
from dataclasses import dataclass
from importlib import import_module


@dataclass(frozen=True)
class AnalysisMethod:
    module: str
    reader: str
    ensemble_input: bool = False

    def analyze(self, settings):
        module = import_module('.' + self.module, package=__package__)
        result = getattr(module, self.reader)(settings if self.ensemble_input else settings['input'])
        return module, result


ANALYSIS_METHODS = {
    '4-points': AnalysisMethod('external_reorganization', 'read_external_four_points'),
    'vertical-gradient': AnalysisMethod('vertical_gradient', 'read_vertical_gradient'),
    'linear-response': AnalysisMethod('gap_reorganization', 'analyze_gap_ensembles', ensemble_input=True),
}
