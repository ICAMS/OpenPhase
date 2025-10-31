$BUILD_INFO_HEADER = "../include/BuildInfo.h"
$GIT_COMMIT_SHA = git rev-parse HEAD
$BUILD_TIME = Get-Date -Format "yyyy-MM-dd HH:mm:ss"

Write-Host "Creating build info header"

# Ensure the directory exists
$directory = Split-Path -Path $BUILD_INFO_HEADER -Parent
if (!(Test-Path -Path $directory)) {
    New-Item -ItemType Directory -Path $directory -Force | Out-Null
}

# Create and write to the file
$headerContent = @"
#ifndef BUILD_INFO_H
#define BUILD_INFO_H
#define GIT_COMMIT_SHA "$GIT_COMMIT_SHA"
#define BUILD_TIME "$BUILD_TIME"
#endif
"@

$headerContent | Set-Content -Path $BUILD_INFO_HEADER

Write-Host "---------------------------------"
