"""cbiokit: helpers for cBioPortal exports.

alterations  driver (OncoKB) events and matrices from 'alterations across samples' exports

The public API is re-exported at the top level::

    import cbiokit as cbk
    ex = cbk.read_alteration_export("alterations_across_samples.tsv")
"""

from .alterations import (AlterationExport, alteration_matrix, concat_exports, driver_event_matrix,
                          driver_gene_matrix, fusion_events, fusion_matrix, fusion_summary, panel_genes_from_export,
                          read_alteration_export)

__version__ = "0.5.1"
