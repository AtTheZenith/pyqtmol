# -*- mode: python ; coding: utf-8 -*-

a = Analysis(
    ['src\\main.py'],
    pathex=['.', 'c:\\Users\\AtTheZenith\\Documents\\Github\\pyqtmol\\.venv\\Lib\\site-packages'],
    binaries=[],
    datas=[],
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='pyqtmol',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    icon=['src\\resources\\icon.ico'],
)

# This is the new section that handles the "onedir" and "contents_directory"
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name='pyqtmol', # This is the main folder name in /dist
    contents_directory='bin' # <--- This is the subfolder with all the packages
)