# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import copy_metadata

datas=[('images','images'),('style.qss','.')]
datas += copy_metadata('scanpy', recursive=True)


block_cipher = None


a = Analysis(
    ['HomeScreen.py','DifferentialGeneAnalysis.py','MainWindowFunctions.py','ThreadHandling.py','VisualizationPopup.py'],
    pathex=[],
    binaries=[],
    datas=datas,
    hiddenimports=['scanpy','opencv-python','PyQt5','leidenalg','openpyxl','anndata'],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='CellVisualization',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='CellVisualization',
)
