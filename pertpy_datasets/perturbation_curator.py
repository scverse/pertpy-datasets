from typing import TYPE_CHECKING, Literal

import anndata as ad

if TYPE_CHECKING:
    from lnschema_core import Record


class _PerturbationValidatorUnavailable:
    """Curator for perturbation data."""

    def __init__(self):
        raise RuntimeError(
            "PerturbationValidator can only be instantiated if connected to a lamindb instance."
        )


# Nested try because django might not be installed
try:
    from django.core.exceptions import ImproperlyConfigured

    try:
        import bionty as bt
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
                    verbosity: The verbosity level.
                    cxg_schema_version: The CELLxGENE schema version to curate against.
                    using_key: A reference LaminDB instance.
                """
                PT_DEFAULT_VALUES = CellxGeneFields.OBS_FIELD_DEFAULTS | {
                    "cell_line": "unknown",
                    "genetic_perturbations": "",
                    "compound_perturbations": "",
                    "environmental_perturbations": "",
                    "combination_perturbations": "",
                }

                PT_CATEGORICALS = CellxGeneFields.OBS_FIELDS | {
                    "cell_line": bt.CellLine.name,
                    "genetic_perturbations": wl.GeneticPerturbation.name,
                    "compound_perturbations": wl.CompoundPerturbation.name,
                    "environmental_perturbations": wl.EnvironmentalPerturbation.name,
                    "combination_perturbations": wl.CombinationPerturbation.name,
                }

                PT_SOURCES: dict[str, Record] = {
                    # "depmap_id": bt.Source.using(using_key)
                    # .filter(name="depmap")
                    # .first(),
                    "cell_line": bt.Source.using(using_key)
                    .filter(name="depmap")
                    .first(),
                    # "compound": bt.Source.using(using_key)
                    # .filter(entity="wetlab.Compound", name="chebi")
                    # .first(),
                }

                self.organism = organism

                # Set the Compound source to chebi; we don't want output if the source has already been set
                with logger.mute():
                    chebi_source = bt.Source.filter(
                        entity="wetlab.Compound", name="chebi"
                    ).first()
                    if not chebi_source:
                        wl.Compound.add_source(
                            bt.Source.filter(entity="Drug", name="chebi").first()
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

            def validate(self) -> bool:
                """Validates the AnnData object against cellxgene and pertpy's requirements."""
                return super().validate()

    except ImproperlyConfigured:
        PerturbationCurator = _PerturbationValidatorUnavailable  # type: ignore

except ImportError:
    PerturbationCurator = _PerturbationValidatorUnavailable  # type: ignore
