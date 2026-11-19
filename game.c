/* ========================================================================
   Strand Open - a genome explorer for anyone's DNA data.

   The unity build: the platform layer (window, settings, processes),
   then the shared header, then each module. See README.md for what it
   does and pipeline/ for how a person's genome database is made.
   ======================================================================== */

#define _GNU_SOURCE   /* strcasestr and friends, also in older C libraries' headers */
#include "src/platform.c"
#include "src/app.h"
#include "src/fonts.c"
#include "src/gui.c"
#include "src/db.c"
#include "src/genome.c"
#include "src/charts.c"
#include "src/geo.c"
#include "src/reads.c"
#include "src/links.c"
#include "src/people.c"
#include "src/page_overview.c"
#include "src/page_browser.c"
#include "src/page_variants.c"
#include "src/page_gene.c"
#include "src/page_health.c"
#include "src/page_pgx.c"
#include "src/page_raw.c"
#include "src/page_ancestry.c"
#include "src/page_rare.c"
#include "src/page_rarity.c"
#include "src/page_brief.c"
#include "src/page_diet.c"
#include "src/report.c"
#include "src/pdf.c"
#include "src/page_report.c"
#include "src/app.c"
