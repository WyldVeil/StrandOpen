/* ========================================================================
   src/win32.c - the Windows side of Strand Open's system calls (see
   win32.h). Compiled on its own for the Windows build only.
   ======================================================================== */
#define WIN32_LEAN_AND_MEAN
#define _WIN32_WINNT 0x0A00
#include <windows.h>
#include <shellapi.h>
#include <shlobj.h>
#include <commdlg.h>
#include <io.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "win32.h"

static void slashes(char *s) { for (; *s; s++) if (*s == '\\') *s = '/'; }
static void backslashes(char *s) { for (; *s; s++) if (*s == '/') *s = '\\'; }

bool w32ExePath(char *out, size_t n) {
    DWORD k = GetModuleFileNameA(NULL, out, (DWORD)n);
    if (k == 0 || k >= n) return false;
    slashes(out);
    return true;
}

bool w32Folder(int which, char *out, size_t n) {
    char p[MAX_PATH];
    int csidl = which == W32_CONFIG ? CSIDL_APPDATA : which == W32_DATA ? CSIDL_LOCAL_APPDATA : CSIDL_PERSONAL;
    if (SHGetFolderPathA(NULL, csidl | CSIDL_FLAG_CREATE, NULL, 0, p) != S_OK) return false;
    snprintf(out, n, "%s", p);
    slashes(out);
    return true;
}

void w32RuntimeEnv(const char *runtime, const char *cacheDir) {
    char rt[MAX_PATH * 2], v[32768];
    snprintf(rt, sizeof rt, "%s", runtime);
    backslashes(rt);
    char *old = getenv("PATH");
    snprintf(v, sizeof v, "%s\\tools;%s\\jre\\bin;%s\\env;%s\\env\\Library\\bin;%s\\env\\Scripts;%s", rt, rt, rt, rt, rt, old ? old : "");
    SetEnvironmentVariableA("PATH", v);
    _putenv_s("PATH", v);
    struct { const char *k; char v[MAX_PATH * 2 + 64]; } e[6];
    snprintf(e[0].v, sizeof e[0].v, "%s\\tools\\plugins", rt);         e[0].k = "BCFTOOLS_PLUGINS";
    snprintf(e[1].v, sizeof e[1].v, "%s\\tools\\busybox.exe", rt);     e[1].k = "STRANDOPEN_BUSYBOX";
    snprintf(e[2].v, sizeof e[2].v, "1");                               e[2].k = "PYTHONNOUSERSITE";
    snprintf(e[3].v, sizeof e[3].v, "1");                               e[3].k = "PYTHONUTF8";
    snprintf(e[4].v, sizeof e[4].v, "1");                               e[4].k = "PYTHONDONTWRITEBYTECODE";
    snprintf(e[5].v, sizeof e[5].v, "%s", cacheDir ? cacheDir : "");   e[5].k = "NUMBA_CACHE_DIR";
    for (int i = 0; i < 6; i++) {
        if (!e[i].v[0]) continue;
        SetEnvironmentVariableA(e[i].k, e[i].v);
        _putenv_s(e[i].k, e[i].v);
    }
}

/* one argument quoted the way the C runtime splits a command line again */
static void quoteArg(const char *a, char *out, size_t n, size_t *len) {
    bool need = !*a || strpbrk(a, " \t\n\v\"") != NULL;
    size_t k = *len;
#define PUT(c) do { if (k + 1 < n) out[k++] = (c); } while (0)
    if (k) PUT(' ');
    if (!need) { for (const char *p = a; *p; p++) PUT(*p); *len = k; out[k] = 0; return; }
    PUT('"');
    for (const char *p = a;; p++) {
        size_t bs = 0;
        while (*p == '\\') { bs++; p++; }
        if (!*p) { for (size_t i = 0; i < bs * 2; i++) PUT('\\'); break; }
        if (*p == '"') { for (size_t i = 0; i < bs * 2 + 1; i++) PUT('\\'); PUT('"'); }
        else { for (size_t i = 0; i < bs; i++) PUT('\\'); PUT(*p); }
    }
    PUT('"');
#undef PUT
    out[k] = 0;
    *len = k;
}

bool w32Spawn(char *const argv[], const char *outFile, int *fdOut, unsigned *pidOut, void **hOut) {
    static char cmd[32768];
    size_t len = 0;
    cmd[0] = 0;
    for (int i = 0; argv[i]; i++) quoteArg(argv[i], cmd, sizeof cmd, &len);
    char exe[MAX_PATH * 2];
    snprintf(exe, sizeof exe, "%s", argv[0]);
    backslashes(exe);

    SECURITY_ATTRIBUTES sa = { sizeof sa, NULL, TRUE };
    HANDLE nul = CreateFileA("NUL", GENERIC_READ | GENERIC_WRITE, FILE_SHARE_READ | FILE_SHARE_WRITE, &sa, OPEN_EXISTING, 0, NULL);
    HANDLE out = INVALID_HANDLE_VALUE, rd = INVALID_HANDLE_VALUE;
    if (outFile) {
        out = CreateFileA(outFile, GENERIC_WRITE, FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE, &sa, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
    } else if (!CreatePipe(&rd, &out, &sa, 0)) {
        out = INVALID_HANDLE_VALUE;
    } else {
        SetHandleInformation(rd, HANDLE_FLAG_INHERIT, 0);
    }
    if (out == INVALID_HANDLE_VALUE) { if (nul != INVALID_HANDLE_VALUE) CloseHandle(nul); return false; }

    STARTUPINFOEXA six;
    memset(&six, 0, sizeof six);
    six.StartupInfo.cb = sizeof six;
    six.StartupInfo.dwFlags = STARTF_USESTDHANDLES;
    six.StartupInfo.hStdInput = nul;
    six.StartupInfo.hStdOutput = out;
    six.StartupInfo.hStdError = outFile ? out : nul;
    /* the child inherits these handles and no others: otherwise it would also get
       every other inheritable handle open at that moment (another reader's pipe,
       this program's own output) and hold them open after we close ours */
    HANDLE list[2] = { nul, out };
    int nl = nul != INVALID_HANDLE_VALUE ? 2 : 1;
    if (nl == 1) list[0] = out;
    SIZE_T asz = 0;
    InitializeProcThreadAttributeList(NULL, 1, 0, &asz);
    LPPROC_THREAD_ATTRIBUTE_LIST al = (LPPROC_THREAD_ATTRIBUTE_LIST)HeapAlloc(GetProcessHeap(), 0, asz);
    DWORD ext = 0;
    if (al && InitializeProcThreadAttributeList(al, 1, 0, &asz)) {
        if (UpdateProcThreadAttribute(al, 0, PROC_THREAD_ATTRIBUTE_HANDLE_LIST, list, nl * sizeof(HANDLE), NULL, NULL)) {
            six.lpAttributeList = al;
            ext = EXTENDED_STARTUPINFO_PRESENT;
        }
    }
    PROCESS_INFORMATION pi;
    /* no console window for it (or for the console programs it starts); left
       on its own when it is the import, so closing the window does not end it */
    DWORD flags = CREATE_NO_WINDOW | ext | (outFile ? CREATE_NEW_PROCESS_GROUP | CREATE_BREAKAWAY_FROM_JOB : 0);
    BOOL ok = CreateProcessA(exe, cmd, NULL, NULL, TRUE, flags, NULL, NULL, &six.StartupInfo, &pi);
    if (!ok && outFile) ok = CreateProcessA(exe, cmd, NULL, NULL, TRUE, flags & ~(DWORD)CREATE_BREAKAWAY_FROM_JOB, NULL, NULL, &six.StartupInfo, &pi);
    if (six.lpAttributeList) DeleteProcThreadAttributeList(al);
    if (al) HeapFree(GetProcessHeap(), 0, al);
    CloseHandle(out);
    if (nul != INVALID_HANDLE_VALUE) CloseHandle(nul);
    if (!ok) { if (rd != INVALID_HANDLE_VALUE) CloseHandle(rd); return false; }
    CloseHandle(pi.hThread);
    if (pidOut) *pidOut = pi.dwProcessId;
    if (outFile) {
        CloseHandle(pi.hProcess);
        if (hOut) *hOut = NULL;
        if (fdOut) *fdOut = -1;
        return true;
    }
    if (hOut) *hOut = pi.hProcess; else CloseHandle(pi.hProcess);
    if (fdOut) *fdOut = _open_osfhandle((intptr_t)rd, _O_RDONLY | _O_BINARY);
    return true;
}

int w32Wait(void *h, bool stop) {
    if (!h) return -1;
    if (stop) TerminateProcess((HANDLE)h, 1);
    WaitForSingleObject((HANDLE)h, INFINITE);
    DWORD code = (DWORD)-1;
    if (!GetExitCodeProcess((HANDLE)h, &code)) code = (DWORD)-1;
    CloseHandle((HANDLE)h);
    return (int)code;
}

bool w32Alive(unsigned pid) {
    HANDLE h = OpenProcess(PROCESS_QUERY_LIMITED_INFORMATION, FALSE, pid);
    if (!h) return false;
    DWORD code = 0;
    bool alive = GetExitCodeProcess(h, &code) && code == STILL_ACTIVE;
    CloseHandle(h);
    return alive;
}

double w32FreeBytes(const char *path) {
    ULARGE_INTEGER avail;
    return GetDiskFreeSpaceExA(path, &avail, NULL, NULL) ? (double)avail.QuadPart : 0;
}

bool w32Replace(const char *from, const char *to) {
    for (int i = 0; i < 20; i++) {     /* a reader may hold the file for a moment */
        if (MoveFileExA(from, to, MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH)) return true;
        Sleep(50);
    }
    return false;
}

bool w32PickFile(void *owner, char *out, size_t n) {
    char buf[MAX_PATH * 4] = "";
    OPENFILENAMEA o;
    memset(&o, 0, sizeof o);
    o.lStructSize = sizeof o;
    o.hwndOwner = (HWND)owner;
    o.lpstrFilter = "DNA files\0*.txt;*.csv;*.tsv;*.zip;*.vcf;*.vcf.gz;*.gz\0All files\0*.*\0";
    o.lpstrFile = buf;
    o.nMaxFile = sizeof buf;
    o.lpstrTitle = "Choose a DNA file";
    o.Flags = OFN_EXPLORER | OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST | OFN_NOCHANGEDIR;
    if (!GetOpenFileNameA(&o)) return false;
    snprintf(out, n, "%s", buf);
    slashes(out);
    return true;
}

bool w32PickFolder(void *owner, const char *title, char *out, size_t n) {
    CoInitializeEx(NULL, COINIT_APARTMENTTHREADED);
    BROWSEINFOA b;
    memset(&b, 0, sizeof b);
    char name[MAX_PATH];
    b.hwndOwner = (HWND)owner;
    b.pszDisplayName = name;
    b.lpszTitle = title;
    b.ulFlags = BIF_RETURNONLYFSDIRS | BIF_NEWDIALOGSTYLE;
    PIDLIST_ABSOLUTE id = SHBrowseForFolderA(&b);
    if (!id) return false;
    char p[MAX_PATH];
    bool ok = SHGetPathFromIDListA(id, p);
    CoTaskMemFree(id);
    if (!ok) return false;
    snprintf(out, n, "%s", p);
    slashes(out);
    return true;
}

void w32Open(const char *pathOrUrl) {
    char p[MAX_PATH * 4];
    snprintf(p, sizeof p, "%s", pathOrUrl);
    if (!strstr(p, "://")) backslashes(p);
    ShellExecuteA(NULL, "open", p, NULL, NULL, SW_SHOWNORMAL);
}

long long w32Pread(int fd, void *buf, size_t n, long long off) {
    HANDLE h = (HANDLE)_get_osfhandle(fd);
    if (h == INVALID_HANDLE_VALUE) return -1;
    OVERLAPPED o;
    memset(&o, 0, sizeof o);
    o.Offset = (DWORD)(off & 0xffffffffLL);
    o.OffsetHigh = (DWORD)(off >> 32);
    DWORD got = 0;
    if (!ReadFile(h, buf, (DWORD)n, &got, &o)) return GetLastError() == ERROR_HANDLE_EOF ? 0 : -1;
    return (long long)got;
}

bool w32ToolsCanRead(const char *path) {
    wchar_t w[MAX_PATH * 4];
    if (!MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, path, -1, w, MAX_PATH * 4)) return false;
    char cp[16] = "";
    UINT acp = GetLocaleInfoA(LOCALE_SYSTEM_DEFAULT, LOCALE_IDEFAULTANSICODEPAGE, cp, sizeof cp) ? (UINT)atoi(cp) : 1252;
    if (acp == 0 || acp == CP_UTF8) return true;
    char out[MAX_PATH * 8];
    BOOL lost = FALSE;
    if (!WideCharToMultiByte(acp, WC_NO_BEST_FIT_CHARS, w, -1, out, sizeof out, NULL, &lost)) return false;
    return !lost;
}

bool w32Relaunch(char *const argv[]) {
    static char cmd[32768];
    size_t len = 0;
    for (int i = 0; argv[i]; i++) quoteArg(argv[i], cmd, sizeof cmd, &len);
    char exe[MAX_PATH * 2];
    snprintf(exe, sizeof exe, "%s", argv[0]);
    backslashes(exe);
    STARTUPINFOA si;
    memset(&si, 0, sizeof si);
    si.cb = sizeof si;
    PROCESS_INFORMATION pi;
    if (!CreateProcessA(exe, cmd, NULL, NULL, FALSE, 0, NULL, NULL, &si, &pi)) return false;
    CloseHandle(pi.hThread);
    CloseHandle(pi.hProcess);
    return true;
}
