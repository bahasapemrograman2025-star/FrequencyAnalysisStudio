# Rainfall Frequency Analysis Studio — Build & Usage Guide

This package contains the full source code for a desktop application
(Windows-friendly, but also runs on macOS/Linux from source) that performs
rainfall / annual-discharge frequency analysis using 6 statistical
distributions, and produces fully formatted native Excel reports with
charts.

Version 3.0 adds:
- **Two data sources, combinable:** import a multi-sheet Excel workbook
  and/or type data manually inside the app (multiple series, just like
  adding sheets in Excel — but you never have to leave the app).
- **A unified Series Manager** that lists every series you've loaded,
  no matter where it came from.
- **A redesigned, dashboard-style UI** with a sidebar, cards, colored
  status badges and a step-by-step workflow (Data Source → Series
  Manager → Process & Results).

---

## 1. What's in this folder

| File                     | Purpose                                                      |
|--------------------------|----------------------------------------------------------------|
| `rainfall_app.py`        | The complete application (analysis engine + GUI).            |
| `requirements.txt`       | Python packages needed to run / build the app.                |
| `build_exe.bat`          | One-click Windows script that builds a standalone `.exe`.    |
| `make_sample_input.py`   | Optional helper that generates `sample_input.xlsx` for testing. |
| `app_icon.ico`           | *(optional, not included)* — drop your own icon here with this exact name to brand the `.exe`. |

---

## 2. Two ways to use this project

### Option A — Run it directly with Python (fastest, works on any OS)

This does **not** produce a `.exe`, but lets you use the application
immediately, on Windows, macOS or Linux.

1. Install **Python 3.10 or newer** from <https://www.python.org/downloads/>.
   - On Windows, tick **"Add Python to PATH"** during installation.
   - `tkinter` (the GUI toolkit) ships with the official Windows and
     macOS installers automatically — no extra step needed there.
     On Linux, if you get a `tkinter not found` error, install it with
     your package manager, e.g. `sudo apt install python3-tk`.
2. Open a terminal / Command Prompt in this folder.
3. Install the dependencies:
   ```
   pip install -r requirements.txt
   ```
4. Run the app:
   ```
   python rainfall_app.py
   ```

### Option B — Build a standalone Windows `.exe` (no Python needed to run it afterwards)

**Important:** `.exe` files are Windows-native binaries. PyInstaller
does **not** cross-compile, so this step **must be run on a Windows
machine** (or a Windows virtual machine). If you are on macOS/Linux,
either use a Windows PC/VM for this step, or use Option A instead.

1. Copy this whole folder onto a Windows computer.
2. Make sure Python 3.10+ is installed on that Windows machine (see
   Option A, step 1).
3. **Double-click `build_exe.bat`** (or right-click → "Run"). Do **not**
   try to run it as `python build_exe.bat` — that is incorrect syntax;
   it is a batch script, not a Python script. Just run it directly.
4. The script will automatically:
   - Create an isolated virtual environment (`build_env`)
   - Install all required packages, including PyInstaller
   - Package everything into one file with PyInstaller
5. When it finishes, your application will be at:
   ```
   dist\RainfallFrequencyAnalysis.exe
   ```
6. That single `.exe` file is fully standalone — copy it anywhere
   (a USB drive, another PC, a shared network folder) and double-click
   to run it. The target machine does **not** need Python installed.

**Optional — custom icon:** place an `app_icon.ico` file in this same
folder before running `build_exe.bat` and the script will automatically
use it as the `.exe` icon. If it's not present, the build still works,
just with the default icon.

**Optional — cleaning up:** after a successful build, you can delete the
`build_env`, `build`, and `dist\*.spec`-related temp folders if you want
to save disk space; the only file you actually need to keep and share
is `dist\RainfallFrequencyAnalysis.exe`.

---

## 3. Troubleshooting the build

| Symptom | Likely cause / fix |
|---|---|
| `'python' is not recognized...` | Python isn't on PATH. Reinstall Python and tick "Add Python to PATH", or use the full path to `python.exe`. |
| Build finishes but the `.exe` won't start / closes instantly | Try building without `--windowed` temporarily (edit `build_exe.bat`, remove `--windowed`) so a console window shows the error message. |
| Windows Defender / antivirus flags the `.exe` | This is a common false-positive with PyInstaller `--onefile` builds (the file self-extracts at runtime, which looks similar to some malware behavior). You can add an exception, or switch to `--onedir` instead of `--onefile` in `build_exe.bat` for a folder-based build that antivirus tools tend to trust more. |
| First launch is slow | `--onefile` builds unpack themselves into a temp folder every time they start. This is normal; later launches on the same session are faster. Use `--onedir` if startup speed matters more than having a single file. |
| `ModuleNotFoundError` when building | Delete the `build_env` folder and re-run `build_exe.bat` to get a clean virtual environment. |

---

## 4. How to use the application

### Step 1 — Data Source
- **Import from Excel:** click "Browse Excel File..." and pick a
  workbook where each sheet has two columns named exactly `Year` and
  `Rainfall (mm)`. Every sheet becomes one series automatically. See
  `sample_input.xlsx` (generate it with `python make_sample_input.py`)
  for a working example.
- **Manual Entry:** click "+ New Manual Series", give it a name, and a
  spreadsheet-style editor opens where you can:
  - Add rows one at a time (`+ Add Row`)
  - Edit or delete a selected row
  - **Paste/Bulk Add** many rows at once (one `year, value` pair per
    line — commas, tabs, semicolons or plain spaces are all accepted)
  - See live statistics (N, mean, std dev, Cv, Cs) update as you type

You can mix both sources freely — imported and manually typed series
all live together in the same workspace.

### Step 2 — Series Manager
A single table listing every series currently loaded, its source,
row count, quick statistics and readiness status. From here you can
rename, edit, or remove any series (or clear everything and start
over).

### Step 3 — Process & Results
- Choose an output folder.
- Click "▶ Process All Series" — every *ready* series (5+ data rows)
  is analyzed using all 6 distributions, with full data-quality
  testing, goodness-of-fit testing (Chi-Square & Kolmogorov–Smirnov),
  and native Excel charts.
- Each series produces its own `.xlsx` report, plus one combined
  `Best_Return_Period_Summary.xlsx` recommending the best-fit
  distribution per series.
- Watch progress and the live log; open the result folder directly
  from the app when done.

---

## 5. Notes

- The analytical engine (distribution fitting, goodness-of-fit tests,
  data-quality tests, and the Excel report layout) is unchanged from
  the original tool — only the data-input workflow and the interface
  were redesigned.
- All on-screen text, log messages and generated reports are in
  English.
