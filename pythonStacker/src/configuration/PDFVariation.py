from src.configuration.Uncertainty import Uncertainty


def generate_pdfvariations(dict_entry):
    ret = [PDFVariation("PDF", dict_entry, i) for i in range(100)]
    return ret


class PDFVariation(Uncertainty):
    def __init__(self, name, dict_entry, instance):
        super().__init__(name, dict_entry)
        self._isFlat = False
        self.type = "weight"
        self._correlated_process = True
        self.correlated_years = True

        self.name = f"PDF_{instance}"
        self.pretty_name = f"PDF_{instance}"
        self.technical_name = f"pdf_{instance}"
        self.weight_key_up = f"PDF_{instance}"
        self.weight_alias_up = f"pdfVariations[:, {instance}]"

        # TODO: implement interpretation of this None
        self.weight_key_down = None
