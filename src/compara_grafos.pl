#!/usr/bin/perl -w

#
# compara_grafos.pl - Programa para verificacao da diferenca simetrica entre uma PGN original e a sua arquitetura recuperada
#
# Argumentos:
# 
# -i=arquivo_rede.txt : arquivo que contem uma descricao de uma PGN
# -o=arquivo_estimada.txt : arquivo que contem uma lista de melhores preditores estimados
# -p=n : calcular o grafo estimado considerando n preditores para cada gene-alvo (default:3)
# -v : modo verbose (mostra em detalhes o grafo original e o estimado)
#
# Devolve, na saida padrao, uma tabela texto, contendo a lista dos melhores preditores para cada gene,
# ordenados em ordem descrescente de pontuacao.
# 
# Copyleft M.S.Reis, 13 de outubro de 2009
#

use strict;

# Variaveis
#
my $verbose = 0;
my $num_preditores = 3;
my $nome_arq_rede = "";
my $nome_arq_rede_estimada = "";

# Le os argumentos
#
foreach(@ARGV){
  if ($_ =~ /\-v/){
    $verbose = 1;
  }
  elsif ($_ =~ /\-i\=(.*)/){
    $nome_arq_rede = $1;
  }
  elsif ($_ =~ /\-o\=(.*)/){
    $nome_arq_rede_estimada = $1;
  }
  elsif ($_ =~ /\-p\=(\d+)/){
    $num_preditores = $1;
    if ($1 <= 0){
      print "Erro! Numero de preditores deve ser inteiro positivo!\n";
      erroParametro();
    }
  }
  else{
    erroParametro();
  }
}

# Abre o arquivo de arquitetura (real) da rede (PGN)
#
if (($nome_arq_rede eq "") or (!open(ARQ1,$nome_arq_rede))){
  print "\nNome de arquivo ($nome_arq_rede) invalido ou arquivo inexistente!\n";
  erroParametro();
}

# Carrega o arquivo de rede
#
my %preditor;    # hash com todos os preditores

while(<ARQ1>){
  chomp $_;
  if (!($_ =~ /^\#.*/)){   # ignora comentarios   
    my @lista = split (" ", $_);
    my $chave = shift @lista;
    shift @lista;  # joga fora 
    my @lista2;
    foreach(@lista){
      $_ =~ /\S*\d+\.(\S+)/;
      push(@lista2, $1);
    }
    $preditor{$chave} = [ @lista2 ];
  }
}
close(ARQ1);

# Abre o arquivo de arquitetura (estimada) da rede (PGN)
#
if (($nome_arq_rede_estimada eq "") or (!open(ARQ2,$nome_arq_rede_estimada))){
  print "\nNome de arquivo ($nome_arq_rede_estimada) invalido ou arquivo inexistente!\n";
  erroParametro();
}

my %preditor_estimado;
my $gene = "";
my $count = 1;

while(<ARQ2>){
  chomp $_;
  if ($_ =~ /^Predicoes\s+de\s+(\w+).*/){ 
    $gene = $1;
    $count = 1;
  }
  elsif(($_ =~ /^(\w+)\s+\:\s+\S+/) and ($count <= $num_preditores)){
    push (@{$preditor_estimado{$gene}}, $1);
    $count++;
  }
}
close(ARQ2);


if($verbose){
  foreach my $gene (sort keys %preditor){
    print "$gene: \tReal ->\t";
    foreach (sort @{$preditor{$gene}}){
      print "$_ ";
    }
    print " \t\tPredito -> ";
    foreach (sort @{$preditor_estimado{$gene}}){
      print "$_ ";
    }
    print "\n";
  }
}

$count = 0;

foreach my $gene (sort keys %preditor){
  my %lista_pred;
  foreach (sort @{$preditor{$gene}}){
    $lista_pred{$_} = 0;
  }
  foreach (sort @{$preditor_estimado{$gene}}){
    if(defined($lista_pred{$_})){
      $lista_pred{$_} = 1;
    }
    else{
      $count++;
    }
  }
  foreach my $chave (keys %lista_pred){
    if($lista_pred{$chave} == 0){
      $count++;
    }
  }
}

if($verbose){
  print "\nDiferenca simetrica entre os arcos da rede original e da estimada: $count\n\n";
}
else{
  print "$count\n";
}

# Fim do "main"
#
exit 0;




#            #
# Subrotinas #
#            #

# Subrotina que exibe a sintaxe do programa e retorna ao sistema.
#
sub erroParametro{

  die "\n\nSintaxe: ./compara_grafos.pl [-v | -p=n] -i=rede.txt -o=rede_estimada.txt\n\n";

}


