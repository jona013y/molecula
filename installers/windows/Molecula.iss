#define MyAppName "Molecula"
#define MyAppVersion "1.0.13"
#define MyAppPublisher "Molecula"

[Setup]
AppId={{8DA58436-357A-46EA-984A-1197D0BBA4A7}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher={#MyAppPublisher}
DefaultDirName={localappdata}\Programs\Molecula
DefaultGroupName={#MyAppName}
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64
OutputDir=..\..\dist\windows
OutputBaseFilename=Molecula-Windows-Setup
SetupIconFile=..\..\assets\Molecula.ico
UninstallDisplayIcon={app}\Molecula.exe
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
CloseApplications=yes
RestartApplications=no
Uninstallable=yes

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "Create a desktop shortcut"; GroupDescription: "Additional shortcuts:"; Flags: unchecked

[Files]
Source: "..\..\target\release\molecula.exe"; DestDir: "{app}"; DestName: "Molecula.exe"; Flags: ignoreversion
Source: "..\..\dist\windows\mesa\*.dll"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\..\dist\windows\mesa-notices\*"; DestDir: "{app}\licenses\mesa"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{autoprograms}\Molecula"; Filename: "{app}\Molecula.exe"; WorkingDir: "{app}"
Name: "{autodesktop}\Molecula"; Filename: "{app}\Molecula.exe"; WorkingDir: "{app}"; Tasks: desktopicon

[Run]
Filename: "{app}\Molecula.exe"; Description: "Launch Molecula"; Flags: postinstall nowait skipifsilent

[UninstallDelete]
Type: filesandordirs; Name: "{localappdata}\Molecula"
