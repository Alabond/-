@echo off
cd /d "%~dp0"
chcp 65001 >nul
set "PATH=D:\GAP\runtime\bin;%PATH%"
set "GAP_EXE=D:\GAP\runtime\opt\gap-4.15.1\gap.exe"
"%GAP_EXE%" -q -A -b G2_Pipeline.g
set "RC=%ERRORLEVEL%"
powershell -NoProfile -ExecutionPolicy Bypass -Command "$utf8NoBom = [System.Text.UTF8Encoding]::new($false); $utf8Bom = [System.Text.UTF8Encoding]::new($true); $files = @('../g2(2)实验运行结果.txt'); foreach($f in $files){ if(Test-Path $f){ $txt = [System.IO.File]::ReadAllText($f,$utf8NoBom); [System.IO.File]::WriteAllText($f,$txt,$utf8Bom) } }"
exit /b %RC%
