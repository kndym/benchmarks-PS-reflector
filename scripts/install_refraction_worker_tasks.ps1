$ErrorActionPreference = 'Stop'
$taskRoot = Split-Path -Parent $PSScriptRoot
$taskUser = [System.Security.Principal.WindowsIdentity]::GetCurrent().Name
$pythonConsole = (Get-Command python.exe).Source
$pythonWindowless = Join-Path (Split-Path -Parent $pythonConsole) 'pythonw.exe'
if (-not (Test-Path -LiteralPath $pythonWindowless)) { throw 'pythonw.exe is required for hidden scheduled workers' }
$principal = New-ScheduledTaskPrincipal -UserId $taskUser -LogonType Interactive -RunLevel Limited
$settings = New-ScheduledTaskSettingsSet -ExecutionTimeLimit ([TimeSpan]::Zero) -StartWhenAvailable -AllowStartIfOnBatteries -DontStopIfGoingOnBatteries -RestartCount 999 -RestartInterval (New-TimeSpan -Minutes 1) -MultipleInstances IgnoreNew
$trigger = New-ScheduledTaskTrigger -AtLogOn -User $taskUser
$names = @{sweep='Reflector-EntropicSweep'; debiased='Reflector-DebiasedEvaluation'; plots='Reflector-ComparisonPlots'}
foreach ($role in @('sweep','debiased','plots')) {
    $action = New-ScheduledTaskAction -Execute $pythonWindowless -Argument "-B scripts/supervise_refraction_runs.py --role $role" -WorkingDirectory $taskRoot
    Register-ScheduledTask -TaskName $names[$role] -Action $action -Principal $principal -Settings $settings -Trigger $trigger -Description 'Checkpoint-resuming reflector worker independent of Codex. Stops when all cases are complete.' -Force | Out-Null
}
foreach ($role in @('sweep','debiased','plots')) {
    Start-ScheduledTask -TaskName $names[$role]
}
Get-ScheduledTask -TaskName 'Reflector-*' | Select-Object TaskName,State
