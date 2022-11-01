(* ::Package:: *)

$FKKSRoot = NotebookDirectory[];
$SolutionPath = $FKKSRoot<>"../Solutions/";
CharacterReprRule = {\[Mu]->mu,\[Nu]->nu,\[Delta]->delta,\[Chi]->chi,\[Eta]->eta,\[Omega]->omega,}


AllData = Import/@FileNames[All, $SolutionPath];


FlattenAssociation[ass_Association]:=KeyMap[Last[#]&,ResourceFunction["AssociationKeyFlatten"][ass]]
RemoveUnpythonable[ass_Association]:=DeleteCases[Map[If[ListQ[#]||StringQ[#]||NumberQ[#], #]&, FlattenAssociation[ass]],_Symbol]
ToSaveableDataset[ass_Association]:=Table[{Keys[ass][[i]], ass[Keys[ass][[i]]]}, {i,1,Length@ass}]
SolutionToPythonable[ass_Association]:=ToSaveableDataset[RemoveUnpythonable[ass]]


SolutionToPythonable[AllData[[1]]]
Export["test.dat", %]
Import["test.dat"]


SymbolName[\[Chi]]//InputForm
