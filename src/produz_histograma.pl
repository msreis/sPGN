#!/usr/bin/perl -w

#
# produz_histograma.pl - Programa para a producao de um histograma:
#                        eixo das abscissas sao os instantes de tempo, enquanto o das ordenadas
#                        eh a diferenca simetrica media entre os arcos da arquitetura real de uma PGN
#                        e a arquitetura obtida atraves da estimacao feita utilizando uma serie temporal
#                        simulada dessa PGN.
#
#             Argumentos:
# Obrigatorios:
# -i=arq.txt : nome do arquivo de PGN a ser considerado
# Opcionais:
# -v   : modo verbose (default off)
# -t=n : numero maximo de instantes de tempo considerados (min 100, max 1000, default 100)
# -r=n : nivel de ruido considerado (min 0, max 100, default: 1)
# -p=n : numero de preditores considerados na estimacao (min 1, default: 3)
# -c   : estimacao utilizando COD (default: IMM)
# -w=n : numero de repeticoes por instante de tempo (default: min 1, max 100, default 10)
# 
# O programa considera para todo valor < 100 "andar" de 10 em 10, e entre 100 e 1000 "andar de 100 em 100.
#
# Devolve, na saida padrao, uma tabela texto, onde a primeira coluna eh o numero de instantes de tempo considerados
# e a segunda eh o valor da diferenca simetrica media.
# 
# Copyleft M.S.Reis, 13 de outubro de 2009
#

use strict;

# Variaveis e constantes;
#
my $VETOR_ESTADOS_ALEATORIO = 1;  # verdadeiro; use zero para falso

my $verbose = 0;
my $num_preditores = 3;
my $arq_rede = "";
my $t = 100;
my $r = 1;
my $c = " ";
my $p = 3;
my $w = 10;

# Le os argumentos
#
foreach(@ARGV){
  if ($_ =~ /\-v/){
    $verbose = 1;
  }
  elsif ($_ =~ /\-i\=(.*)/){
    $arq_rede = $1;
  }
  elsif ($_ =~ /\-p\=(\d+)/){
    $num_preditores = $1;
    if ($1 <= 0){
      print "Erro! Numero de preditores deve ser inteiro positivo!\n";
      erroParametro();
    }
  }
  elsif ($_ =~ /\-r\=(\d+)/){
    $r = $1;
    if (($1 == 0) or ($1>100)){
      print "Erro! Ruido eh um valor entre [0, 100]!\n";
      erroParametro();
    }
  }
  elsif ($_ =~ /\-w\=(\d+)/){
    $w = $1;
    if (($1 == 0) or ($1>100)){
      print "Erro! Repeticoes eh um valor inteiro entre [0, 100]!\n";
      erroParametro();
    }
  }
  elsif ($_ =~ /\-t\=(\d+)/){
    $t = $1;
    if (($1 < 10) or ($1>1000)){
      print "Erro! Tamanho maximo da serie temporal eh um valor inteiro entre [10, 1000]!\n";
      erroParametro();
    }
  }
  elsif ($_ =~ /\-c/){
    $c = " -c ";
  }
  else{
    erroParametro();
  }
}

if($arq_rede eq ""){
  print "Arquivo de PGN invalido!\n";
  erroParametro();
}

open(OUT, ">$arq_rede.hist") or die "Erro ao abrir o arquivo de saida!\n";

for (my $i = 10; $i <= 100; $i = $i + 10){
  geraSerie($i);
}

for (my $i = 200; $i <= $t; $i = $i + 100){
  geraSerie($i);
}

close(OUT);

# Fim do "main"
#
exit 0;



#            #
# Subrotinas #
#            #

# Subrotina que exibe a sintaxe do programa e retorna ao sistema.
#
sub erroParametro{

  die "\n\nSintaxe: ./produz_histograma.pl [-v] -p=n -r=n -c -t=n -w=n -i=rede.txt\n\n";

}


# subrotina que calcula um E[X] para uma serie temporal com i instantes de tempo, gravando o resultado
# no arquivo de saida
#
sub geraSerie{

  my $i = shift @_;
  
  my $E_X;  # E[X]

  # calcula o E[X] para o i atual, tirando a media com w realizacoes de X
  #
  for (my $j = 1; $j <= $w; $j++){
    
    my $a = " ";

    if($VETOR_ESTADOS_ALEATORIO){
      $a = " -a ";
    }

    # simula a serie temporal
    #
    print "Simulando uma PGN com $i instantes de tempo... ";
    system("./sPGN.pl -t=" . $i . "-r=" . $r . $a . "-i=" . $arq_rede . " > $arq_rede.out.sPGN.tmp");

    # estima a rede utilizando o metodo escolhido (IMM ou CoD)
    #
    print "[done]\nEstimando uma PGN com a serie produzida... ";
    system("./ePGN.pl " . $c . "-i=$arq_rede.out.sPGN.tmp -p=" . $p . " > $arq_rede.out.ePGN.tmp");

    # calcula a diferenca simetrica entre a arquitetura da rede real e a da estimada
    #
    print"[done]\nCalculando a diferenca simetrica e apagando os arquivos temporarios... ";
    system("./compara_grafos.pl -i=" . $arq_rede . " -p=" . $p . " -o=$arq_rede.out.ePGN.tmp > $arq_rede.out.compara_grafos.tmp");

    # soma em E[X] o valor de X atual;
    #
    open(GRAFO, "$arq_rede.out.compara_grafos.tmp") or die "Erro ao abrir a saida do compara_grafos.pl!\n";
    $_ = <GRAFO>;
    $_ =~ /(\d+)/;
    $E_X += $1;
    close(GRAFO);

    # apaga os arquivos temporarios
    #
    system ("rm $arq_rede.out.sPGN.tmp");
    system ("rm $arq_rede.out.ePGN.tmp");
    system ("rm $arq_rede.out.compara_grafos.tmp");

    print "[done]\n";

  }

  # computa o E[X] e escreve no arquivo de saida
  #
  $E_X = $E_X / $w;
  printf OUT "%d %3.3f\n", $i, $E_X;

  # fim da subrotina 
}

