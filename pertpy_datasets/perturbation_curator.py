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


def standardize_value_unit(value: str, is_dose: bool = True) -> tuple[str, str]:
    """Standardize a string containing a value and unit into a tuple of (float, unit).

    Handles various input formats and standardizes unit case.
    Always formats numbers with exactly one decimal place.
    """
    # Handle non-string or empty inputs
    if not isinstance(value, str) or not value.strip():
        return None

    # Strip whitespace
    value = str(value).strip()

    # Regular expression to match number and unit
    match = re.match(r"^(\d*\.?\d{0,1})\s*([a-zA-ZμµΜ]+)$", value)

    if not match:
        raise ValueError(
            f"Invalid format: {value}. Expected format: number with max 1 decimal place + unit (e.g., '10.0nM' or '5.0h')"
        )

    number, unit = match.groups()
    # Format to always have exactly 1 decimal place
    formatted_number = f"{float(number):.1f}"
    # Standardize unit case
    if is_dose:
        standardized_unit = standardize_dose_unit(unit)
    else:
        standardized_unit = unit.lower()
    return formatted_number, standardized_unit


def validate_dose_unit(unit: str) -> bool:
    """Validate if the dose unit is acceptable.

    Only accepts specific molar concentration units with exact case:
    - nM: nanomolar
    - μM: micromolar (with Greek mu)
    - mM: millimolar
    - M: molar
    """
    valid_units = {"nM", "μM", "mM", "M"}
    return unit in valid_units


def validate_time_unit(unit: str) -> bool:
    """Validate if the time unit is acceptable.

    Strictly accepts:
    - s: seconds
    - m: minutes
    - h: hours
    - d: days
    - y: years
    """
    valid_units = {"h", "m", "s", "d", "y"}
    return unit in valid_units


def validate_dose_time(df: pd.DataFrame) -> None:
    """Validate pert_dose and pert_time columns.

    Empty strings and NA values are allowed.
    Other invalid values will raise errors.
    """
    errors = []

    # Validate dose column
    for idx, value in df["pert_dose"].items():
        # Skip NA/empty/'nan' values
        if pd.isna(value):
            continue

        if isinstance(value, (int, float)):
            errors.append(f"Row {idx} - Missing unit for dose: {value}. Must include a unit (nM, μM, mM, M)")
            continue

        try:
            _, unit = standardize_value_unit(value, is_dose=True)
            if not validate_dose_unit(unit):
                errors.append(f"Row {idx} - Invalid dose unit: {unit}. Must be exactly one of: nM, μM, mM, M")
        except ValueError as e:
            errors.append(f"Row {idx} - {str(e)}")

    # Validate time column
    for idx, value in df["pert_time"].items():
        # Skip NA/empty/'nan' values
        if pd.isna(value):
            continue

        if isinstance(value, (int, float)):
            errors.append(f"Row {idx} - Missing unit for time: {value}. Must include a unit (h, m, s, d, y)")
            continue

        try:
            _, unit = standardize_value_unit(value, is_dose=False)
            if not validate_time_unit(unit):
                errors.append(
                    f"Row {idx} - Invalid time unit: {unit}. Use 'h' for hours, 'm' for minutes, 's' for seconds, 'd' for days, or 'y' for years"
                )
        except ValueError as e:
            errors.append(f"Row {idx} - {str(e)}")

    if errors:
        raise ValueError("Validation errors found:\n" + "\n".join(errors))


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
        from lamin_utils import logger
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
                        raise ValueError(
                            "Either 'pert_target' or both 'pert_name' and 'pert_type' must be present in adata.obs"
                        )
                else:
                    if "pert_name" not in adata.obs.columns:
                        logger.warning("No 'pert' column found in adata.obs, will only curate 'pert_target'")
                    elif "pert_type" not in adata.obs.columns:
                        raise ValueError("Both 'pert' and 'pert_type' must be present in adata.obs")

                # check for invalid perturbation types
                if "pert_type" in adata.obs.columns:
                    data_pert_types = set(adata.obs["pert_type"].unique())
                    invalid_pert_types = data_pert_types - set(PERT_COLUMNS.keys())
                    if invalid_pert_types:
                        raise ValueError(
                            f"invalid pert_type found: {invalid_pert_types}!\n    → allowed values: {set(PERT_COLUMNS.keys())}"
                        )
                    for pert_type in data_pert_types:
                        adata.obs[PERT_COLUMNS[pert_type]] = adata.obs["pert_name"].where(
                            adata.obs["pert_type"] == pert_type, None
                        )

                # validate dose and time columns
                if pert_dose:
                    if "pert_dose" not in adata.obs.columns:
                        raise ValueError(
                            "pert_dose column is not found!    → please provide a pert_dose column or set pert_dose=False"
                        )
                    else:
                        if not ln.Feature.filter(name="pert_dose").exists():
                            ln.Feature(name="pert_dose", dtype="str").save()
                if pert_time:
                    if "pert_time" not in adata.obs.columns:
                        raise ValueError(
                            "pert_time column is not found!    → please provide a pert_time column or set pert_time=False"
                        )
                    else:
                        if not ln.Feature.filter(name="pert_time").exists():
                            ln.Feature(name="pert_time", dtype="str").save()
                validate_dose_time(adata.obs)

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

            def validate(self) -> bool:
                """Validates the AnnData object against cellxgene and pertpy's requirements."""
                return super().validate()

            def standardize_dose_time(self) -> None:
                """Standardize pert_dose and pert_time columns in place.

                NA values are preserved as is.

                Args:
                    df: Input DataFrame with 'pert_dose' and 'pert_time' columns

                Returns:
                    DataFrame with standardized values in pert_dose and pert_time columns
                """
                standardized_df = self._adata.obs.copy()

                # Process dose column
                for idx, value in self._adata.obs["pert_dose"].items():
                    if pd.isna(value) or (isinstance(value, str) and (not value.strip() or value.lower() == "nan")):
                        standardized_df.at[idx, "pert_dose"] = None
                        continue
                    try:
                        num, unit = standardize_value_unit(value, is_dose=True)
                        if validate_dose_unit(unit):
                            standardized_df.at[idx, "pert_dose"] = f"{num}{unit}"
                    except ValueError:
                        continue

                # Process time column
                for idx, value in self._adata.obs["pert_time"].items():
                    if pd.isna(value) or (isinstance(value, str) and (not value.strip() or value.lower() == "nan")):
                        standardized_df.at[idx, "pert_time"] = None
                        continue
                    try:
                        num, unit = standardize_value_unit(value, is_dose=False)
                        if validate_time_unit(unit):
                            standardized_df.at[idx, "pert_time"] = f"{num}{unit}"
                    except ValueError:
                        continue

                self._adata.obs = standardized_df

    except ImproperlyConfigured:
        PerturbationCurator = _PerturbationValidatorUnavailable  # type: ignore

except ImportError:
    PerturbationCurator = _PerturbationValidatorUnavailable  # type: ignore
