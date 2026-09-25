This directory contains stub files for testing `get_gtdb.pl` without
having to handle very large data.

`get_gtdb.pl download|all --test` copies these files instead of downloading
(selected by download label, from this directory next to the script, so any
working directory works):

| Download label      | Stub                                  |
|---------------------|---------------------------------------|
| `bac_markers`       | `bac_markers.tar.gz`                  |
| `arc_markers`       | `arc_markers.tar.gz`                  |
| `bac_taxonomy`      | `bac_taxonomy.tsv.gz`                 |
| `arc_taxonomy`      | `arc_taxonomy.tsv.gz`                 |
| `bac_metadata`      | `bac_metadata.tsv.gz`                 |
| `arc_metadata`      | `arc_metadata.tsv.gz`                 |
| `tk_database_parts` | `gtdbtk_dummy.tar.gz.part_{aa,ab,ac}` |
| `tk_database`       | the three parts concatenated          |

The marker archives use the GTDB layout `<archive>/{faa,fna}/<marker>` (two
markers each, none shared between bacteria and archaea); the tables hold one
genome each (neither is a GTDB representative); the GTDB-Tk parts form a
`.tar.gz` of `gtdbtk_dummy/{aa,ab,ac}/*.txt`.

`t/get_gtdb.t` runs the script end to end on these stubs (and on archives it
builds itself, e.g. for markers shared between bacteria and archaea):

    prove -I. t/get_gtdb.t
