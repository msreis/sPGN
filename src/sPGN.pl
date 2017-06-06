#!/usr/bin/perl -w

#
# sPGN.pl - Programa para simulacao da dinamica de uma Probabilistic Genetic Network (PGN).
#
# Argumentos:
# 
# -i=arquivo.txt : arquivo que contem uma especificacao de arquitetura de uma PGN
# -t=n : simular a execucao de n instantes de tempo (min:1, max:1000, default: 100)
# -r=m : simular a execucao com m porcento de ruido (min:0, max:100, default: 1)
# -a   : iniciar a simulacao com um vetor de estados aleatorio.
# --nome_gene=t : da um pulso no gene "nome_gene", no instante de tempo t 
#
# Devolve, na saida padrao, uma tabela texto, contendo os valores dos genes ao longo do periodo simulado.
# 

# Modificacoes:
# 11.out.09 - modificada a simulacao para uma PGN "standard" (paper da malaria), removendo os thresholds.

# Formato de um arquivo de entrada:
#
# primeira coluna: nome do gene-alvo i  
# segunda coluna: valor inicial do gene-alvo ( S_i[0] )
# terceira coluna em diante: genes preditores de x, acompanhados de seu respectivo coeficiente ( a_i_j )
#
# Exemplo:  x00 1  -2.x01 3.x02 1.x03  
#           Significa que x00 e' o gene-alvo, inicia com valor "1", e e' predito por x01 (com coeficiente -2),
#           por x02 (com coeficiente 3) e por x03 (com coeficiente 1).
# 

# Copyleft M.S.Reis, 14 de agosto 2009
#

use strict;

# Definicoes
#
my $MAX_NUM_ITER = 1000;
my $MIN_NUM_ITER = 1;
my $MAX_NOISE = 100;
my $MIN_NOISE = 0;
my $MAX_VALUE = 2;   # discretizado em 3 estados: 0, 1 ou 2
my $MIN_VALUE = 0;

my %pulso;

# Argumentos
#
my $t = 100;
my $r = 1;
my $nome_arq = "";
my $vetor_aleatorio = 0;

# Le os argumentos
#
foreach(@ARGV){
  if ($_ =~ /\-t\=(\d+)/){
    if(($1 < $MIN_NUM_ITER) or ($1 > $MAX_NUM_ITER)){
      erroParametro();
    } 
    else{
      $t = $1;
    }
  }
  elsif ($_ =~ /\-r\=(\d+)/){
    if(($1 < $MIN_NOISE) or ($1 > $MAX_NOISE)){
      erroParametro();
    } 
    else{
      $r = $1;
    }
  }
  elsif ($_ =~ /\-i\=(.*)/){
    $nome_arq = $1;
  }
  elsif ($_ =~ /\-a/){
    $vetor_aleatorio = 1;
  }
  elsif ($_ =~ /\-\-(\S+)\=(\d+)/){
    $pulso{$1} = $2;
  }

  else{
    erroParametro();
  }
}

# Le o arquivo de arquitetura da rede (PGN)
#
if (($nome_arq eq "") or (!open(ARQ,$nome_arq))){
  print "\nNome de arquivo ($nome_arq) invalido ou arquivo inexistente!\n";
  erroParametro();
}

# Variaveis:
#
my %S;    # hash com os estados dos genes no instante de tempo t
my %a;    # hash com a lista de preditores de cada gene

while(<ARQ>){
  chomp $_;
  if (!($_ =~ /^\#.*/)){   # ignora comentarios    
    my @lista = split (" ", $_);
    my $gene = shift @lista;
    my $S_0 = shift @lista;
    $a{$gene} = [ @lista ];  # guardo a lista de preditores em um hash, com o alvo como chave
    if ($vetor_aleatorio){
      $S{$gene}->[0] = int(rand($MAX_VALUE + 1));
    }
    else{
      $S{$gene}->[0] = $S_0;
    }
  }
}
close(ARQ);

# Simulando a execucao de t instantes de tempo (t=0 exclusive), com r porcento de ruido:
#
for(my $i=1; $i<=$t; $i++){

  foreach my $gene (sort keys %S){

    my $d = 0;

    foreach (@{$a{$gene}}){
      $_ =~  /(\S+)\.(\S+)/;
      $d += $1 * $S{$2}->[$i-1];    #  pred. do gene i: SUM_j ( a_ij * S_j(t-1) ),
    }                               #  onde a_ij  eh "1" ou "-1" para todo preditor j do gene i
    
    # calculo do Y:
    #
    my $y = $S{$gene}->[$i-1] + $d;

    if ($y > 2){
      $y = 2;
    }
    elsif ($y < 0){
      $y = 0;
    }

    # calculo dos ruidos
    #
    my $ruido = rand(100);

    if($ruido >= $r){
      $S{$gene}->[$i] = $y;
    }
    else{   # com r porcento de chances, cai em um estado diferente de y:
      my $k;
      do{
	$k = int(rand($MAX_VALUE + 1));
	# print "d=$y, k=$k\n";
      }while($k == $y);
      $S{$gene}->[$i] = $k;
    }
    
    # se $gene tem um pulso no instante $i, atualiza S{gene}->[i] com o valor do pulso
    if ((defined($pulso{$gene})) and ($i == $pulso{$gene})){
      $S{$gene}->[$i] = 1;
    }
    elsif ((defined($pulso{$gene})) and ($i-1 == $pulso{$gene})){
      $S{$gene}->[$i] = 0;
    }
  }
}

# Imprime os resultados na saida padrao, um gene por linha
#
foreach my $gene (sort keys %S){
  print "$gene ";
  for (my $i=0; $i<=$t; $i++){
    printf "%d ", $S{$gene}->[$i];
  }
  print "\n";
}

exit 0;
#
# Fim do "main"
#


#            #
# Subrotinas #
#            #

# Subrotina que exibe a sintaxe do programa e retorna ao sistema.
#
sub erroParametro{

  die "\n\nSintaxe: ./sPGN.pl [-t=n] [-r=m] [-a] -i=arquitetura.txt\n\n";

}
