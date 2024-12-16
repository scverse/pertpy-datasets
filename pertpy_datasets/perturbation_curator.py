import re
from typing import TYPE_CHECKING, Literal

import anndata as ad
import pandas as pd

if TYPE_CHECKING:
    from lnschema_core import Record


def standardize_dose_unit(unit: str) -> str:
    """Standardize dose unit to standard formats."""
    unit_map = {
        "nm": "nM",
        "NM": "nM",
        "um": "μM",
        "UM": "μM",
        "μm": "μM",
        "μM": "μM",
        "µm": "μM",
        "µM": "μM",
        "mm": "mM",
        "MM": "mM",
        "m": "M",
        "M": "M",
    }
    return unit_map.get(unit, unit)


def validate_dose_unit(unit: str) -> bool:
    """Validate if the dose unit is acceptable.

    Only accepts specific molar concentration units with exact case:
    - nM: nanomolar
    - μM: micromolar (with Greek mu)
    - mM: millimolar
    - M: molar.
    """
    valid_units = {
        "nM",
        "μM",
        "µM",  # Alternative for μM
        "mM",
        "M",
    }
    return unit in valid_units


def validate_time_unit(unit: str) -> bool:
    """Validate if the time unit is acceptable.

    Strictly accepts ONLY:
    - s: seconds
    - m: minutes
    - h: hours
    - d: days
    - y: years
    """
    valid_units = {"h", "m", "s", "d", "y"}
    return unit == unit.lower() and unit in valid_units


def parse_value_unit(value: str, is_dose: bool = True) -> tuple[str, str]:
    """Parse a string containing a value and unit into a tuple without standardization."""
    if not isinstance(value, str) or not value.strip():
        return None

    value = str(value).strip()
    match = re.match(r"^(\d*\.?\d{0,1})\s*([a-zA-ZμµΜ]+)$", value)

    if not match:
        raise ValueError(
            f"Invalid format: {value}. Expected format: number with max 1 decimal place + unit (e.g., '10.0nM' or '5.0h')"
        )

    number, unit = match.groups()
    formatted_number = f"{float(number):.1f}"

    if is_dose and not validate_dose_unit(unit):
        raise ValueError(f"Invalid dose unit: {unit}. Must be exactly one of: nM, μM, mM, M (case sensitive)")
    elif not is_dose and not validate_time_unit(unit):
        raise ValueError(f"Invalid time unit: {unit}. Must be exactly one of: h, m, s, d, y (single letter only)")

    return formatted_number, unit


def standardize_value_unit(value: str, is_dose: bool = True) -> tuple[str, str]:
    """Standardize a string containing a value and unit."""
    if not isinstance(value, str) or not value.strip():
        return None

    value = str(value).strip()
    match = re.match(r"^(\d*\.?\d{0,1})\s*([a-zA-ZμµΜ]+)$", value)

    if not match:
        raise ValueError(
            f"Invalid format: {value}. Expected format: number with max 1 decimal place + unit (e.g., '10.0nM' or '5.0h')"
        )

    number, unit = match.groups()
    formatted_number = f"{float(number):.1f}"

    if is_dose:
        standardized_unit = standardize_dose_unit(unit)
        if not validate_dose_unit(standardized_unit):
            raise ValueError(f"Invalid dose unit: {unit}. Must be convertible to one of: nM, μM, mM, M")
    else:
        # For time, strip any extra characters (hr -> h, min -> m, etc)
        if unit.startswith("hr"):
            standardized_unit = "h"
        elif unit.startswith("min"):
            standardized_unit = "m"
        elif unit.startswith("sec"):
            standardized_unit = "s"
        else:
            standardized_unit = unit[0].lower()

        if not validate_time_unit(standardized_unit):
            raise ValueError(f"Invalid time unit: {unit}. Must be convertible to one of: h, m, s, d, y")

    return formatted_number, standardized_unit


def validate_dose(values: pd.Series) -> list:
    """Validate pert_dose values with strict case checking."""
    errors = []

    for idx, value in values.items():
        if pd.isna(value):
            continue

        if isinstance(value, (int, float)):
            errors.append(f"Row {idx} - Missing unit for dose: {value}. Must include a unit (nM, μM, mM, M)")
            continue

        try:
            parse_value_unit(value, is_dose=True)
        except ValueError as e:
            errors.append(f"Row {idx} - {str(e)}")

    return errors


def validate_time(values: pd.Series) -> list:
    """Validate pert_time values."""
    errors = []

    for idx, value in values.items():
        if pd.isna(value):
            continue

        if isinstance(value, (int, float)):
            errors.append(f"Row {idx} - Missing unit for time: {value}. Must include a unit (h, m, s, d, y)")
            continue

        try:
            parse_value_unit(value, is_dose=False)
        except ValueError as e:
            errors.append(f"Row {idx} - {str(e)}")

    return errors


class _PerturbationValidatorUnavailable:
    """Curator for perturbation data."""

    def __init__(self):
        raise RuntimeError("PerturbationValidator can only be instantiated if connected to a lamindb instance.")


# Nested try because django might not be installed
try:
    from django.core.exceptions import ImproperlyConfigured

    try:
        import bionty as bt
        import lamindb as ln
        import wetlab as wl
        from cellxgene_lamin import CellxGeneFields, Curate
        from lamin_utils import colors, logger
        from lamindb.core.exceptions import ValidationError
        from lamindb_setup.core.types import UPathStr
        from lnschema_core.types import FieldAttr

        class PerturbationCurator(Curate):
            """Curator flow for Perturbation data."""

            def __init__(
                self,
                adata: ad.AnnData | UPathStr,
                var_index: FieldAttr = bt.Gene.ensembl_gene_id,
                organism: Literal["human", "mouse"] = "human",
                pert_dose: bool = True,
                pert_time: bool = True,
                *,
                verbosity: str = "hint",
                cxg_schema_version: Literal["5.0.0", "5.1.0"] = "5.1.0",
                using_key: str = "laminlabs/pertpy-datasets",
            ):
                """Curator flow for Perturbation data.

                Args:
                    adata: Path to or AnnData object to curate against the CELLxGENE schema.
                    var_index: The registry field for mapping the ``.var`` index.
                    organism: The organism name. CELLxGENE restricts it to 'human' and 'mouse' and therefore so do we.
                    pert_dose: Whether to curate the 'pert_dose' column.
                    pert_time: Whether to curate the 'pert_time' column.
                    verbosity: The verbosity level.
                    cxg_schema_version: The CELLxGENE schema version to curate against.
                    using_key: A reference LaminDB instance.
                """
                PERT_COLUMNS = {
                    "compound": "pert_compound",
                    "genetic": "pert_genetic",
                }
                PT_DEFAULT_VALUES = CellxGeneFields.OBS_FIELD_DEFAULTS | {
                    "cell_line": "unknown",
                    "pert_target": "unknown",
                    "pert_genetic": "",
                    "pert_compound": "",
                    # "compound_perturbation": "",
                    # "environmental_perturbation": "",
                    # "combination_perturbation": "",
                }

                PT_CATEGORICALS = CellxGeneFields.OBS_FIELDS | {
                    "cell_line": bt.CellLine.name,
                    "pert_target": wl.PerturbationTarget.name,
                    "pert_genetic": wl.GeneticPerturbation.name,
                    "pert_compound": wl.Compound.name,
                    # "compound_perturbation": wl.CompoundPerturbation.name,
                    # "environmental_perturbation": wl.EnvironmentalPerturbation.name,
                    # "combination_perturbation": wl.CombinationPerturbation.name,
                }

                PT_SOURCES: dict[str, Record] = {
                    "cell_line": bt.Source.using(using_key).filter(name="depmap").first(),
                    "pert_compound": bt.Source.using(using_key).filter(entity="wetlab.Compound", name="chebi").first(),
                }

                self.organism = organism

                # Set the Compound source to chebi; we don't want output if the source has already been set
                with logger.mute():
                    chebi_source = bt.Source.filter(entity="wetlab.Compound", name="chebi").first()
                    if not chebi_source:
                        wl.Compound.add_source(bt.Source.filter(entity="Drug", name="chebi").first())

                # check required perturbation columns
                if "pert_target" not in adata.obs.columns:
                    if "pert_name" not in adata.obs.columns or "pert_type" not in adata.obs.columns:
                        raise ValidationError(
                            "either 'pert_target' or both 'pert_name' and 'pert_type' must be present in adata.obs"
                        )
                else:
                    if "pert_name" not in adata.obs.columns:
                        logger.warning("no 'pert' column found in adata.obs, will only curate 'pert_target'")
                    elif "pert_type" not in adata.obs.columns:
                        raise ValidationError("both 'pert' and 'pert_type' must be present in adata.obs")

                # check for invalid perturbation types
                if "pert_type" in adata.obs.columns:
                    data_pert_types = set(adata.obs["pert_type"].unique())
                    invalid_pert_types = data_pert_types - set(PERT_COLUMNS.keys())
                    if invalid_pert_types:
                        raise ValidationError(
                            f"invalid pert_type found: {invalid_pert_types}!\n    → allowed values: {set(PERT_COLUMNS.keys())}"
                        )
                    for pert_type in data_pert_types:
                        adata.obs[PERT_COLUMNS[pert_type]] = adata.obs["pert_name"].where(
                            adata.obs["pert_type"] == pert_type, None
                        )

                # validate dose and time columns
                if pert_dose:
                    if "pert_dose" not in adata.obs.columns:
                        raise ValidationError(
                            "pert_dose column is not found!    → please provide a pert_dose column or set pert_dose=False"
                        )
                if pert_time:
                    if "pert_time" not in adata.obs.columns:
                        raise ValidationError(
                            "pert_time column is not found!    → please provide a pert_time column or set pert_time=False"
                        )

                super().__init__(
                    adata=adata,
                    var_index=var_index,
                    categoricals=PT_CATEGORICALS,
                    using_key=using_key,
                    defaults=PT_DEFAULT_VALUES,
                    verbosity=verbosity,
                    organism=self.organism,
                    extra_sources=PT_SOURCES,
                    schema_version=cxg_schema_version,
                )

                self._pert_time = pert_time
                self._pert_dose = pert_dose

            def validate(self) -> bool:
                """Validates the AnnData object against cellxgene and pertpy's requirements."""
                validated = super().validate()

                # validate dose and time columns
                if self._pert_dose:
                    if not ln.Feature.filter(name="pert_dose").exists():
                        ln.Feature(name="pert_dose", dtype="str").save()
                    dose_errors = validate_dose(self._adata.obs["pert_dose"])
                    if dose_errors:
                        errors_print = "\n    ".join(dose_errors)
                        logger.warning(
                            f"invalid pert_dose values found!\n    {errors_print}\n    → run {colors.cyan('standardize_dose_time()')}"
                        )
                        validated &= False
                if self._pert_time:
                    if not ln.Feature.filter(name="pert_time").exists():
                        ln.Feature(name="pert_time", dtype="str").save()
                    time_errors = validate_time(self._adata.obs["pert_time"])
                    if time_errors:
                        errors_print = "\n    ".join(time_errors)
                        logger.warning(
                            f"invalid pert_time values found!\n    {errors_print}\n   → run {colors.cyan('standardize_dose_time()')}"
                        )
                        validated &= False
                self._validated = validated
                return validated

            def standardize_dose_time(self) -> pd.DataFrame:
                """Standardize pert_dose and pert_time columns.

                This function will attempt to standardize units to their correct case:
                - Dose units will be standardized to nM, μM, mM, M
                - Time units will be standardized to lowercase h, m, s, d, y

                NA values are preserved as is.

                Args:
                    df: DataFrame with 'pert_dose' and 'pert_time' columns

                Returns:
                    DataFrame with standardized values
                """
                standardized_df = self._adata.obs.copy()

                # Process dose column
                if "pert_dose" in self._adata.obs.columns:
                    for idx, value in self._adata.obs["pert_dose"].items():
                        if pd.isna(value) or (isinstance(value, str) and (not value.strip() or value.lower() == "nan")):
                            standardized_df.at[idx, "pert_dose"] = None
                            continue
                        try:
                            num, unit = standardize_value_unit(value, is_dose=True)
                            standardized_df.at[idx, "pert_dose"] = f"{num}{unit}"
                        except ValueError:
                            continue

                # Process time column
                if "pert_time" in self._adata.obs.columns:
                    for idx, value in self._adata.obs["pert_time"].items():
                        if pd.isna(value) or (isinstance(value, str) and (not value.strip() or value.lower() == "nan")):
                            standardized_df.at[idx, "pert_time"] = None
                            continue
                        try:
                            num, unit = standardize_value_unit(value, is_dose=False)
                            standardized_df.at[idx, "pert_time"] = f"{num}{unit}"
                        except ValueError:
                            continue

                self._adata.obs = standardized_df

    except ImproperlyConfigured:
        PerturbationCurator = _PerturbationValidatorUnavailable  # type: ignore

except ImportError:
    PerturbationCurator = _PerturbationValidatorUnavailable  # type: ignore
