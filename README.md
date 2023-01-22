![](https://github.com/joaosds/3body_rel/tree/main/img/Triquette.gif)

### Base - Three Body


- [ ] Arrumar git, documentação geral e README;
- [ ] Arrumar escalas;
- [ ] Generalizar funções para n corpos de maneira geral;
- [ ] Consistentemente entender a tolerância na integração e como o erro influencia o resultado final.
- [ ]  Arrumar método de integração, com conservação de energia (integração simplética) - julia; Leap-frog não leva em conta uma força dependente da velocidade, por exemplo.
- [ ] Solução para quando planetas se chocam. Por que na modificação do Felipe isso já é contemplado? Verificar se uma das partículas está emitindo NAN e não sendo registrada. Ideia: Epsilon dependente do raio;
- [ ] Gráfico do Anel de partículas sob influência das GW
- [ ] Perfil das GW possui 3 fases. Nosso código captura bem a primeira. A segunda necessita relatividade numérica (talvez usar einsteinpy), e a terceira pode ser calculada com teoria de perturbação. 

- [ ] Existe um banco de dados de estrelas binárias para a 2a parte. Se tivessemos para 3 corpos, talvez um ML? 

### Base - Black Hole Image 
- [ ] Arrumar região do meio para foto de buraco negro (elipse);

### Ideias futuras

- [ ] Não lembro da ideia inicial envolvendo astrofísica 
- [ ] Reconstrução no diagrama de fases (dado parâmetros iniciais) da evolução do sistema caótico. Compreensão dos objetos antes da colisão.
- [ ] Seção de choque de buracos negros (?)

### Fase Final
- [ ] Pensar ao fim como melhorar o design do pyqt;
- [ ] Artigo;
- [ ] Arrumar sites 

### Próxima reunião geral (28/01/23)

- [ ] Talvez traduzir o codigo pra Julia


# Links
1. Old website https://mcgill3body.github.io/
2. Old rep https://github.com/mcgill3body/mcgill3body.github.io/tree/master/sourcecode

![](https://github.com/joaosds/3body_rel/tree/main/img/logo_project_name_site.png)
