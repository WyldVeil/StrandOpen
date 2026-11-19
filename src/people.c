/* ========================================================================
   src/people.c - people (profiles), importing a person's file, switching.

   Data root (setting so.root, or STRANDOPEN_ROOT, or ~/.local/share/
   StrandOpen):
     shared.kv            where the shared reference data lives (ref=, kg=)
     profiles/<slug>/     one person: profile.kv (name, sex, input.1),
                          canon/, work/, db/strand.db, run/progress.kv
   The import is pipeline/run.py, started detached so it carries on if the
   window closes; its progress file is polled here. Opening a person
   restarts the program with --profile <slug>, so nothing from the previous
   person can linger in any page's cache.
   ======================================================================== */

#include <dirent.h>

static char g_profileSlug[64];

typedef struct {
    char slug[64], name[160], sex[16], input[1024], kind[64], state[24], title[120], error[300], skipped[600];
    int done, of;
    bool hasDb;
} Person;
static Person g_people[64];
static int g_peopleN;
static double g_peopleAt;

/* key=value file into a callback-free little table */
static bool kvRead(const char *path, const char *key, char *out, size_t n) {
    FILE *f = fopen(path, "r");
    if (!f) return false;
    char line[2048];
    size_t kl = strlen(key);
    bool ok = false;
    while (fgets(line, sizeof line, f)) {
        if (!strncmp(line, key, kl) && line[kl] == '=') {
            line[strcspn(line, "\r\n")] = 0;
            snprintf(out, n, "%s", line + kl + 1);
            ok = true;
            break;
        }
    }
    fclose(f);
    return ok;
}

static void soResolveRoot(void) {
    const char *e = getenv("STRANDOPEN_ROOT");
    const char *c = cfgGet("so.root", "");
    if (e && *e) snprintf(g_root, sizeof g_root, "%s", e);
    else if (*c) snprintf(g_root, sizeof g_root, "%s", c);
    else {
#ifdef _WIN32
        char base[PATH_MAX];
        if (w32Folder(W32_DATA, base, sizeof base)) snprintf(g_root, sizeof g_root, "%s/StrandOpen", base);
        else snprintf(g_root, sizeof g_root, "%s/data", g_appDir);
#else
        const char *x = getenv("XDG_DATA_HOME"), *h = getenv("HOME");
        if (x && *x) snprintf(g_root, sizeof g_root, "%s/StrandOpen", x);
        else snprintf(g_root, sizeof g_root, "%s/.local/share/StrandOpen", h && *h ? h : "/tmp");
#endif
    }
    char p[1100];
    snprintf(p, sizeof p, "%s/profiles", g_root);
    mkdirs(p, 0700);
    snprintf(p, sizeof p, "%s/shared.kv", g_root);
    if (!kvRead(p, "ref", g_refDir, sizeof g_refDir)) snprintf(g_refDir, sizeof g_refDir, "%s/shared/ref", g_root);
}

static void peopleScan(void) {
    g_peopleN = 0;
    char dir[1100];
    snprintf(dir, sizeof dir, "%s/profiles", g_root);
    DIR *d = opendir(dir);
    if (!d) return;
    struct dirent *e;
    while ((e = readdir(d)) && g_peopleN < (int)NELEMS(g_people)) {
        if (e->d_name[0] == '.') continue;
        Person *p = &g_people[g_peopleN];
        memset(p, 0, sizeof *p);
        char f[1400];
        snprintf(f, sizeof f, "%s/%s/profile.kv", dir, e->d_name);
        if (!kvRead(f, "name", p->name, sizeof p->name)) continue;
        snprintf(p->slug, sizeof p->slug, "%s", e->d_name);
        kvRead(f, "sex", p->sex, sizeof p->sex);
        kvRead(f, "input.1", p->input, sizeof p->input);
        snprintf(f, sizeof f, "%s/%s/canon/ingest.kv", dir, e->d_name);
        char fmt[32] = "", vendor[64] = "", cov[32] = "";
        kvRead(f, "format", fmt, sizeof fmt);
        kvRead(f, "vendor", vendor, sizeof vendor);
        kvRead(f, "coverage", cov, sizeof cov);
        snprintf(p->kind, sizeof p->kind, "%s", !strcmp(fmt, "chip") ? sstr("%s chip", *vendor ? vendor : "Genotyping") :
                                                  !strcmp(fmt, "gvcf") ? "Whole genome (gVCF)" :
                                                  !strcmp(fmt, "vcf") ? (!strcmp(cov, "genome") ? "Whole genome (VCF)" : "Variant calls (VCF)") : "");
        snprintf(f, sizeof f, "%s/%s/run/progress.kv", dir, e->d_name);
        char v[32] = "";
        kvRead(f, "state", p->state, sizeof p->state);
        kvRead(f, "title", p->title, sizeof p->title);
        kvRead(f, "error", p->error, sizeof p->error);
        kvRead(f, "skipped", p->skipped, sizeof p->skipped);
        if (kvRead(f, "done", v, sizeof v)) p->done = atoi(v);
        if (kvRead(f, "of", v, sizeof v)) p->of = atoi(v);
        /* a "running" import whose process is gone was interrupted */
        if (!strcmp(p->state, "running") && kvRead(f, "pid", v, sizeof v) && !procAlive(atol(v)))
            snprintf(p->state, sizeof p->state, "interrupted");
        snprintf(f, sizeof f, "%s/%s/db/strand.db", dir, e->d_name);
        p->hasDb = access(f, R_OK) == 0;
        g_peopleN++;
    }
    closedir(d);
    for (int i = 1; i < g_peopleN; i++)
        for (int j = i; j > 0 && strcasecmp(g_people[j - 1].name, g_people[j].name) > 0; j--) {
            Person t = g_people[j]; g_people[j] = g_people[j - 1]; g_people[j - 1] = t;
        }
    g_peopleAt = GetTime();
}

/* a pipeline script, detached: a grandchild in its own session, so it survives
   the window; its output goes to logf */
static bool spawnDetached(char *const argv[], const char *logf) {
#ifdef _WIN32
    unsigned pid = 0;
    return w32Spawn(argv, logf, NULL, &pid, NULL);
#else
    pid_t pid = fork();
    if (pid < 0) return false;
    if (pid == 0) {
        setsid();
        if (fork() != 0) _exit(0);
        int fd = open(logf, O_WRONLY | O_CREAT | O_TRUNC, 0600);
        int nul = open("/dev/null", O_RDONLY);
        if (nul >= 0) dup2(nul, 0);
        if (fd >= 0) { dup2(fd, 1); dup2(fd, 2); }
        execv(argv[0], argv);
        _exit(127);
    }
    int st;
    waitpid(pid, &st, 0);
    return true;
#endif
}

/* the bundled Python (runtime/bin) first, then the system's */
static bool pipelinePaths(const char *script, char *py, size_t pyN, char *path, size_t pathN) {
    if (!findTool("python3", py, pyN)) return false;
    snprintf(path, pathN, "%s/pipeline/%s", g_appDir, script);
    if (access(path, R_OK) != 0) snprintf(path, pathN, "%s/../pipeline/%s", g_appDir, script);
    return access(path, R_OK) == 0;
}

static bool importStart(const char *slug) {
    char py[1100] = "", runpy[1100], pdir[1100], logf[1200];
    if (!pipelinePaths("run.py", py, sizeof py, runpy, sizeof runpy)) return false;
    snprintf(pdir, sizeof pdir, "%s/profiles/%s", g_root, slug);
    snprintf(logf, sizeof logf, "%s/run", pdir);
    mkdirs(logf, 0700);
    snprintf(logf, sizeof logf, "%s/run/launcher.txt", pdir);
    /* a fresh status file, so the page does not show the last run's result */
    char pk[1200];
    snprintf(pk, sizeof pk, "%s/run/progress.kv", pdir);
    FILE *f = fopen(pk, "w");
    if (f) { fprintf(f, "state=running\ntitle=Starting\ndone=0\nof=20\n"); fclose(f); }
    char *argv[] = { py, "-B", runpy, "--profile", pdir, NULL };
    return spawnDetached(argv, logf);
}

/* ---- setup: the shared reference data, downloaded once ---- */
typedef struct {
    char state[24], phase[160], item[160], error[300];
    double done, total;
    int step, steps;
    bool essentials, population;
} SetupInfo;
static SetupInfo g_setup;

static void setupScan(void) {
    memset(&g_setup, 0, sizeof g_setup);
    char p[1200], v[64], kgDir[1024];
    snprintf(p, sizeof p, "%s/shared.kv", g_root);
    if (!kvRead(p, "ref", g_refDir, sizeof g_refDir)) snprintf(g_refDir, sizeof g_refDir, "%s/shared/ref", g_root);
    if (!kvRead(p, "kg", kgDir, sizeof kgDir)) snprintf(kgDir, sizeof kgDir, "%s/shared/1kg", g_root);
    char a[1200], b[1200];
    snprintf(a, sizeof a, "%s/GRCh38.fa.fai", g_refDir);
    snprintf(b, sizeof b, "%s/raw/clinvar.vcf.gz", g_refDir);
    g_setup.essentials = access(a, R_OK) == 0 && access(b, R_OK) == 0;
    snprintf(a, sizeof a, "%s/pgen/chr22.pgen", kgDir);
    g_setup.population = access(a, R_OK) == 0;
    snprintf(p, sizeof p, "%s/shared/setup/progress.kv", g_root);
    kvRead(p, "state", g_setup.state, sizeof g_setup.state);
    kvRead(p, "phase", g_setup.phase, sizeof g_setup.phase);
    kvRead(p, "item", g_setup.item, sizeof g_setup.item);
    kvRead(p, "error", g_setup.error, sizeof g_setup.error);
    if (kvRead(p, "done_bytes", v, sizeof v)) g_setup.done = atof(v);
    if (kvRead(p, "total_bytes", v, sizeof v)) g_setup.total = atof(v);
    if (kvRead(p, "step", v, sizeof v)) g_setup.step = atoi(v);
    if (kvRead(p, "steps", v, sizeof v)) g_setup.steps = atoi(v);
    if (!strcmp(g_setup.state, "running") && kvRead(p, "pid", v, sizeof v) && !procAlive(atol(v)))
        snprintf(g_setup.state, sizeof g_setup.state, "interrupted");
}

static bool setupStart(bool population) {
    char py[1100] = "", sp[1100], dir[1200], logf[1300];
    if (!pipelinePaths("setup.py", py, sizeof py, sp, sizeof sp)) return false;
    snprintf(dir, sizeof dir, "%s/shared/setup", g_root);
    mkdirs(dir, 0700);
    snprintf(logf, sizeof logf, "%s/launcher.txt", dir);
    char pk[1300];
    snprintf(pk, sizeof pk, "%s/progress.kv", dir);
    FILE *f = fopen(pk, "w");
    if (f) { fprintf(f, "state=running\nphase=Starting\n"); fclose(f); }
    cfgSet("so.setup_pop", population ? "1" : "0");
    cfgSave();
    char *argv[] = { py, "-B", sp, "--root", g_root, population ? "--population" : NULL, NULL };
    return spawnDetached(argv, logf);
}

/* the card at the top of the People page; returns its height */
static float setupCard(float x, float y, float W) {
    bool running = !strcmp(g_setup.state, "running"), failed = !strcmp(g_setup.state, "failed") || !strcmp(g_setup.state, "interrupted");
    if (g_setup.essentials && g_setup.population && !running) return 0;
    double freeGb = freeBytes(g_root) / 1e9;
    static const char essText[] =
        "Strand Open compares a person's DNA with public reference data, downloaded once and shared by everyone you add. The essentials "
        "(about 3 GB to download, 7 GB on disk) cover health, traits, medicines and lineages. The population reference (about 33 GB more "
        "to download, 12 GB kept) adds ancestry, polygenic scores, rarity and Neanderthal DNA. Nothing about anyone's genome is sent.";
    float essH = txtParaH(F_REG, TH->fsSmall, W - 40, 19, essText, 6);
    float h = running ? 118 : failed ? 150 : !g_setup.essentials ? 46 + essH + 72 : 126;
    Rectangle c = { x, y, W, h };
    gxCard(c);
    drRound((Rectangle){ c.x, c.y + 14, 4, h - 28 }, 2, running ? TH->accent : failed ? TH->bad : TH->warn);
    if (running) {
        txt(F_SEMI, TH->fsH3, c.x + 20, c.y + 16, TH->text, "Setting up the reference data");
        bool dl = g_setup.total > 0 && !strncmp(g_setup.phase, "Downloading", 11);
        txtIn(F_REG, TH->fsSmall, (Rectangle){ c.x + 20, c.y + 46, W - 40, 18 }, TH->dim, 0,
              dl ? sstr("Downloading %s (file %d of %d)", g_setup.item, g_setup.step, g_setup.steps) : *g_setup.phase ? g_setup.phase : "Starting");
        float fr = dl ? fminf(1.0f, (float)(g_setup.done / fmax(1, g_setup.total))) : (float)fmod(GetTime() * 0.25, 1.0);
        gxMeter((Rectangle){ c.x + 20, c.y + 74, W * 0.6f, 8 }, fr, TH->accent);
        if (dl) txt(F_REG, TH->fsTiny + 1, c.x + 30 + W * 0.6f, c.y + 69, TH->faint, sstr("%.1f of %.1f GB", g_setup.done / 1e9, g_setup.total / 1e9));
        txtIn(F_REG, TH->fsTiny + 1, (Rectangle){ c.x + 20, c.y + 92, W - 40, 16 }, TH->faint, 0,
              "It carries on in the background if you close the window, and resumes where it stopped if interrupted.");
        return h;
    }
    if (failed) {
        txt(F_SEMI, TH->fsH3, c.x + 20, c.y + 16, TH->text, "The reference download stopped");
        txtIn(F_REG, TH->fsSmall, (Rectangle){ c.x + 20, c.y + 46, W - 40, 18 }, TH->bad, 0,
              !strcmp(g_setup.state, "interrupted") ? "It was interrupted (the computer slept or was shut down)." : g_setup.error);
        txtIn(F_REG, TH->fsTiny + 1, (Rectangle){ c.x + 20, c.y + 70, W - 40, 16 }, TH->faint, 0,
              "Trying again keeps everything already downloaded.");
        if (gxButtonEx(gxId("setup retry"), (Rectangle){ c.x + 20, c.y + 100, 150, 32 }, "Try again", BTN_PRIMARY, true, NULL))
            setupStart(!strcmp(cfgGet("so.setup_pop", "0"), "1"));
        return h;
    }
    if (!g_setup.essentials) {
        txt(F_SEMI, TH->fsH3, c.x + 20, c.y + 16, TH->text, "First: the reference data");
        txtPara(F_REG, TH->fsSmall, c.x + 20, c.y + 46, W - 40, 19, TH->dim, essText, 6);
        float by = c.y + 46 + essH + 18;
        if (gxButtonEx(gxId("setup ess"), (Rectangle){ c.x + 20, by, 230, 34 }, "Download the essentials", BTN_NORMAL, freeGb > 8,
                       "About 3 GB; health, traits, medicines and lineages"))
            setupStart(false);
        if (gxButtonEx(gxId("setup all"), (Rectangle){ c.x + 262, by, 290, 34 }, "Essentials and population", BTN_PRIMARY, freeGb > 55,
                       "About 36 GB to download; needs 55 GB free while it builds"))
            setupStart(true);
        txtIn(F_REG, TH->fsTiny + 1, (Rectangle){ c.x + 568, by + 9, W - 588, 16 }, freeGb > 55 ? TH->faint : TH->warn, 0,
              freeGb > 55 ? sstr("%.0f GB free here", freeGb) : sstr("%.0f GB free here: the population reference needs 55 GB while it builds", freeGb));
        return h;
    }
    txt(F_SEMI, TH->fsH3, c.x + 20, c.y + 16, TH->text, "Optional: the population reference");
    txtPara(F_REG, TH->fsSmall, c.x + 20, c.y + 46, W - 300, 19, TH->dim,
            "Adds ancestry, polygenic scores, rarity and Neanderthal DNA: about 33 GB to download and 12 GB kept, an hour or so of building "
            "afterwards. Re-import people once it is done.", 3);
    if (gxButtonEx(gxId("setup pop"), (Rectangle){ c.x + W - 260, c.y + 46, 240, 34 }, "Add the population reference", BTN_PRIMARY, freeGb > 50,
                   "Needs about 50 GB free while it builds"))
        setupStart(true);
    txt(F_REG, TH->fsTiny + 1, c.x + W - 260, c.y + 88, freeGb > 50 ? TH->faint : TH->warn, sstr("%.0f GB free here", freeGb));
    return h;
}

static void profileOpen(const char *slug) {
    cfgSet("so.profile", slug);
    cfgSave();
    char *args[] = { "--profile", (char *)slug, NULL };
    relaunch(args);
}

/* ---- adding a person ---- */
static char g_addName[120], g_addPath[1024], g_addMsg[300];
static int g_addSex;               /* 0 not sure, 1 male, 2 female */
static bool g_addOpen;

static void addPick(void) {
#ifdef _WIN32
    char pick[1100];
    if (w32PickFile(GetWindowHandle(), pick, sizeof pick)) snprintf(g_addPath, sizeof g_addPath, "%s", pick);
    return;
#endif
    char tool[1100];
    char *argv[8];
    if (findTool("kdialog", tool, sizeof tool)) {
        char *a[] = { tool, "--getopenfilename", getenv("HOME") ? getenv("HOME") : "/", "*.txt *.csv *.zip *.vcf *.vcf.gz *.gz|DNA files", NULL };
        memcpy(argv, a, sizeof a);
    } else if (findTool("zenity", tool, sizeof tool)) {
        char *a[] = { tool, "--file-selection", "--title=Choose a DNA file", NULL };
        memcpy(argv, a, sizeof a);
    } else {
        snprintf(g_addMsg, sizeof g_addMsg, "No file chooser found: type or paste the path, or drop the file on the window.");
        return;
    }
    Proc p = spawnReader(argv);
    if (!p.pid) return;
    char buf[1100] = "";
    ssize_t n = read(p.fd, buf, sizeof buf - 1);
    reap(&p, false);
    if (n > 0) {
        buf[n] = 0;
        buf[strcspn(buf, "\r\n")] = 0;
        if (*buf) snprintf(g_addPath, sizeof g_addPath, "%s", buf);
    }
}

static void addCreate(void) {
    struct stat st;
    if (!g_setup.essentials) { snprintf(g_addMsg, sizeof g_addMsg, "Download the reference data first (above): the import needs it."); return; }
    if (!g_addName[0]) { snprintf(g_addMsg, sizeof g_addMsg, "Give the person a name (any label you like)."); return; }
    if (stat(g_addPath, &st) != 0 || !S_ISREG(st.st_mode)) { snprintf(g_addMsg, sizeof g_addMsg, "That file cannot be read."); return; }
    if (!pathToolSafe(g_addPath)) {
        snprintf(g_addMsg, sizeof g_addMsg, "The import's tools cannot read this file's folder name on this computer (it has letters outside "
                 "the system's character set). Move the file to a folder with a plain name, such as C:\\DNA, and choose it there.");
        return;
    }
    char slug[16];
    unsigned seed = (unsigned)time(NULL) ^ (unsigned)getpid();
    for (int tries = 0; tries < 50; tries++) {
        for (int i = 0; i < 6; i++) slug[i] = "abcdefghijkmnpqrstuvwxyz23456789"[rand_r(&seed) % 32];
        slug[6] = 0;
        char d[1100];
        snprintf(d, sizeof d, "%s/profiles/%s", g_root, slug);
        if (mkdir(d, 0700) == 0) {
            char pk[1200];
            snprintf(pk, sizeof pk, "%s/profile.kv", d);
            FILE *f = fopen(pk, "w");
            if (!f) break;
            fprintf(f, "name=%s\nsex=%s\ninput.1=%s\n", g_addName, g_addSex == 1 ? "male" : g_addSex == 2 ? "female" : "unknown", g_addPath);
            fclose(f);
            if (!importStart(slug)) snprintf(g_addMsg, sizeof g_addMsg, "The import could not be started (Python missing?).");
            else { g_addOpen = false; g_addName[0] = g_addPath[0] = 0; g_addMsg[0] = 0; }
            peopleScan();
            return;
        }
    }
    snprintf(g_addMsg, sizeof g_addMsg, "Could not create the person's folder under %s.", g_root);
}

/* where everything lives: a folder chooser, then a restart on the new root */
static void rootChange(void) {
    char tool[1100], buf[1100] = "";
    Proc p = { 0 };
#ifdef _WIN32
    (void)tool; (void)p;
    if (!w32PickFolder(GetWindowHandle(), "Choose where Strand Open keeps its data", buf, sizeof buf)) return;
#else
    if (findTool("kdialog", tool, sizeof tool)) {
        char *a[] = { tool, "--getexistingdirectory", g_root, NULL };
        p = spawnReader(a);
    } else if (findTool("zenity", tool, sizeof tool)) {
        char *a[] = { tool, "--file-selection", "--directory", "--title=Choose where Strand Open keeps its data", NULL };
        p = spawnReader(a);
    } else return;
    if (!p.pid) return;
    ssize_t n = read(p.fd, buf, sizeof buf - 1);
    reap(&p, false);
    if (n <= 0) return;
    buf[n] = 0;
    buf[strcspn(buf, "\r\n")] = 0;
#endif
    if (!*buf || !strcmp(buf, g_root) || access(buf, W_OK) != 0) return;
    cfgSet("so.root", sstr("%s/StrandOpen", buf));
    cfgSet("so.profile", "");
    cfgSave();
    relaunch(NULL);
}

static void peopleDrops(void) {
    if (!IsFileDropped()) return;
    FilePathList l = LoadDroppedFiles();
    if (l.count > 0) { snprintf(g_addPath, sizeof g_addPath, "%s", l.paths[0]); g_addOpen = true; }
    UnloadDroppedFiles(l);
}

/* ---- the People page ---- */
static GxScroll g_peopleScroll;
static const char *const g_sexOpts[] = { "Not sure", "Male", "Female" };

static void pagePeople(Rectangle r) {
    if (g_peopleAt == 0 || GetTime() - g_peopleAt > 1.0) { peopleScan(); setupScan(); }
    peopleDrops();
    float top = gxScrollBegin(&g_peopleScroll, r);
    float x = r.x + 32, W = fminf(1100, r.width - 64);
    float y = pageHeader((Rectangle){ r.x, top, r.width, r.height }, "People",
                         "Each person's DNA results, kept apart. Add a person from their raw data file (23andMe, AncestryDNA, MyHeritage, "
                         "FamilyTreeDNA, LivingDNA, or a VCF from whole-genome sequencing); the import runs in the background and keeps going "
                         "if you close the window. Everything stays on this computer.");
    float sh = setupCard(x, y, W);
    if (sh > 0) y += sh + 16;
    for (int i = 0; i < g_peopleN; i++) {
        Person *p = &g_people[i];
        bool running = !strcmp(p->state, "running"), failed = !strcmp(p->state, "failed") || !strcmp(p->state, "interrupted");
        float h = 92 + (failed && *p->error ? 22 : 0);
        Rectangle c = { x, y, W, h };
        gxCard(c);
        bool cur = !strcmp(p->slug, g_profileSlug);
        if (cur) drRound((Rectangle){ c.x, c.y + 14, 4, h - 28 }, 2, TH->accent);
        txt(F_SEMI, TH->fsH3, c.x + 20, c.y + 16, TH->text, p->name);
        const char *sub = sstr("%s%s%s", *p->kind ? p->kind : "Not imported yet", *p->sex && strcmp(p->sex, "unknown") ? " - " : "",
                               *p->sex && strcmp(p->sex, "unknown") ? p->sex : "");
        txtIn(F_REG, TH->fsSmall, (Rectangle){ c.x + 20, c.y + 42, W * 0.5f, 18 }, TH->dim, 0, sub);
        if (running) {
            float fr = p->of > 0 ? (float)p->done / (float)p->of : 0;
            gxMeter((Rectangle){ c.x + 20, c.y + 68, W * 0.45f, 8 }, fr, TH->accent);
            txtIn(F_REG, TH->fsTiny + 1, (Rectangle){ c.x + 30 + W * 0.45f, c.y + 62, W * 0.4f, 18 }, TH->dim, 0,
                  sstr("Importing: %s (%d of %d)", p->title, p->done, p->of));
        } else if (failed) {
            txtIn(F_SEMI, TH->fsSmall, (Rectangle){ c.x + 20, c.y + 64, W - 260, 18 }, TH->bad, 0,
                  !strcmp(p->state, "interrupted") ? "The import was interrupted" : sstr("The import stopped at: %s", p->title));
            if (*p->error) txtIn(F_REG, TH->fsTiny + 1, (Rectangle){ c.x + 20, c.y + 86, W - 260, 16 }, TH->faint, 0, p->error);
        } else if (p->hasDb) {
            if (g_setup.population && strstr(p->skipped, "1000 Genomes"))
                txtIn(F_SEMI, TH->fsTiny + 1, (Rectangle){ c.x + 20, c.y + 66, W - 260, 16 }, TH->warn, 0,
                      "The population reference is ready now: re-import to add ancestry, polygenic scores and rarity.");
            else
                txtIn(F_REG, TH->fsTiny + 1, (Rectangle){ c.x + 20, c.y + 66, W - 260, 16 }, TH->faint, 0,
                      *p->skipped ? sstr("Ready. Not available for this data: %s", p->skipped) : "Ready.");
        }
        float bx = c.x + W - 20;
        if (p->hasDb && !running) {
            bx -= 110;
            if (gxButtonEx(gxIdN(60000 + i), (Rectangle){ bx, c.y + 16, 110, 30 }, cur ? "Showing" : "Open", cur ? BTN_GHOST : BTN_PRIMARY, !cur, NULL))
                profileOpen(p->slug);
        }
        if (!running) {
            bx -= 120;
            if (gxButtonEx(gxIdN(61000 + i), (Rectangle){ bx, c.y + 16, 112, 30 }, failed || !p->hasDb ? "Import again" : "Re-import", BTN_GHOST, g_setup.essentials,
                           "Read the file again and rebuild this person's results"))
                { importStart(p->slug); peopleScan(); }
        }
        y += h + 12;
    }
    if (!g_peopleN) {
        y += txtPara(F_REG, TH->fs, x, y, W, 21, TH->dim, "Nobody here yet. Add the first person below.", 2) + 14;
    }
    /* add a person */
    y += 10;
    if (!g_addOpen) {
        if (gxButtonEx(gxId("people add"), (Rectangle){ x, y, 180, 34 }, "Add a person", BTN_PRIMARY, true, NULL)) g_addOpen = true;
        txt(F_REG, TH->fsSmall, x + 196, y + 9, TH->faint, "or drop a DNA file on this window");
        y += 50;
    } else {
        Rectangle c = { x, y, W, 236 + (g_addMsg[0] ? 24 : 0) };
        gxCard(c);
        txt(F_SEMI, TH->fsH3, c.x + 20, c.y + 16, TH->text, "Add a person");
        txt(F_REG, TH->fsSmall, c.x + 20, c.y + 52, TH->dim, "Name");
        gxTextField(gxId("add name"), (Rectangle){ c.x + 140, c.y + 46, 320, 32 }, g_addName, sizeof g_addName, "Any label, e.g. Mum", false);
        txt(F_REG, TH->fsSmall, c.x + 20, c.y + 96, TH->dim, "Sex");
        gxSegmented(gxId("add sex"), (Rectangle){ c.x + 140, c.y + 90, 320, 30 }, g_sexOpts, 3, &g_addSex);
        txt(F_REG, TH->fsTiny + 1, c.x + 476, c.y + 98, TH->faint, "Not sure: it is read from the data");
        txt(F_REG, TH->fsSmall, c.x + 20, c.y + 140, TH->dim, "DNA file");
        gxTextField(gxId("add path"), (Rectangle){ c.x + 140, c.y + 134, W - 300, 32 }, g_addPath, sizeof g_addPath, "Path to the raw data file (.txt, .csv, .zip, .vcf.gz)", false);
        if (gxButtonEx(gxId("add pick"), (Rectangle){ c.x + W - 150, c.y + 134, 130, 32 }, "Choose...", BTN_NORMAL, true, NULL)) addPick();
        if (g_addMsg[0]) txtIn(F_REG, TH->fsSmall, (Rectangle){ c.x + 140, c.y + 176, W - 160, 18 }, TH->warn, 0, g_addMsg);
        float by = c.y + c.height - 50;
        if (gxButtonEx(gxId("add go"), (Rectangle){ c.x + 140, by, 160, 32 }, "Import", BTN_PRIMARY, true, NULL)) addCreate();
        if (gxButtonEx(gxId("add cancel"), (Rectangle){ c.x + 312, by, 110, 32 }, "Cancel", BTN_GHOST, true, NULL)) { g_addOpen = false; g_addMsg[0] = 0; }
        y += c.height + 16;
    }
    double fb = freeBytes(g_root);
    if (fb > 0)
        txtIn(F_REG, TH->fsTiny + 1, (Rectangle){ x, y + 8, W - 190, 16 }, TH->faint, 0,
              sstr("Data folder: %s (%.0f GB free)", g_root, fb / 1e9));
    if (!pathToolSafe(g_root))
        txtIn(F_REG, TH->fsSmall, (Rectangle){ x, y + 30, W - 190, 40 }, TH->warn, 0,
              "This folder's name has letters the import's tools cannot read on this computer, so imports would fail. "
              "Choose a folder with a plain name (Change folder...), such as C:\\StrandOpen.");
    if (gxButtonEx(gxId("people root"), (Rectangle){ x + W - 170, y, 170, 30 }, "Change folder...", BTN_GHOST, true,
                   "Keep the reference data and everyone's results somewhere else (a bigger disk); the program restarts there"))
        rootChange();
    y += 48;
    gxScrollEnd(&g_peopleScroll, r, y - top);
}
