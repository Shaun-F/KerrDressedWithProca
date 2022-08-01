$FKKSRoot = FileNameJoin[FileNameSplit[ExpandFileName[First[$ScriptCommandLine]]][[1;;-3]]]<>"/";
Import[FileNameJoin[{$FKKSRoot,"Packages", "KerrWithProca.wl"}]];
Import[FileNameJoin[{$FKKSRoot, "Packages", "HelperFunctions.wl"}]];
parameterlist = Import[$FKKSRoot<>"WolframScript/parameters.ini"]//ToExpression;
parameters = KeyMap[ToString[ToExpression[#]]&, parameterlist];

printData = parameters[[{"\[Epsilon]", "\[Mu]Nv", "m", "\[Eta]", "n", "l", "s","\[Chi]", "KMax", "branch"}]];


Print[parameters];
Print[assocToString[printData]];
