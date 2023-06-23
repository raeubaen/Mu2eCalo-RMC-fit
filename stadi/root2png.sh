#!/bin/bash
for var in "$@"
do
export filename=${var%%.*}
root -l $var <<-EOF
TCanvas *c = (TCanvas*)_file0->Get(_file0->GetListOfKeys()->At(0)->GetName());
c->Draw();
c->SaveAs("$filename.png");
EOF

done
