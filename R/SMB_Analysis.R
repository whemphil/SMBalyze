blind.input <- function(path.to.files='./'){
  files=list.files(path = path.to.files,pattern = '[.]tif')
  visible.info=sub('_.*','',files)
  blind.ids=sample(LETTERS[1:length(files)],length(files))
  old.file.ids=files[order(blind.ids)]
  new.file.ids=paste0(blind.ids[order(blind.ids)],'.tif')
  file.visibles=visible.info[order(blind.ids)]
  for (i in 1:length(files)){
    dir.create(path = paste0(path.to.files,LETTERS[i]),recursive = TRUE)
    file.copy(from = paste0(path.to.files,old.file.ids[i]),to = paste0(path.to.files,LETTERS[i],'/',new.file.ids[i]),overwrite = TRUE)
  }
  write.table(x=data.frame('new.file.ids'=new.file.ids,'old.file.ids'=old.file.ids),file = paste0(path.to.files,'file_key.txt'),quote = FALSE,col.names = TRUE,row.names = FALSE,sep = '\t')
  write.table(x=data.frame('new.file.ids'=new.file.ids,'visible.info'=file.visibles),file = paste0(path.to.files,'file_visibles.txt'),quote = FALSE,col.names = TRUE,row.names = FALSE,sep = '\t')
}

id.spots <- function(time.step,path.to.file='./',file.name=NULL,spot.box=6,spot.radius=6,spot.min=NULL,spot.max=Inf,spot.picking='composite'){
  find.spots <- function(data,box.size,low.lim,high.lim,fill.radius,spot.picking){
    lq.calc <- function(input){
      input.2=na.omit(c(input))
      result=(input.2[order(input.2)])[round(0.25*length(input.2))]
      return(result)
    }
    vol.calc <- function(input,x,y,radius){
      lq.calc <- function(input){
        input.2=na.omit(c(input))
        result=(input.2[order(input.2)])[round(0.25*length(input.2))]
        return(result)
      }
      result=sum(input[(x-radius):(x+radius),(y-radius):(y+radius)])-lq.calc(input[(x-radius):(x+radius),(y-radius):(y+radius)])*(2*radius+1)^2
      return(result)
    }
    if (spot.picking=='composite'){
      if (is.null(low.lim)==TRUE){
        low.lim=500
      }
      col.id=t(array(1:dim(data)[2],dim = dim(data)[c(2,1)]))
      row.id=array(1:dim(data)[1],dim = dim(data))
      peeks=array(F,dim = dim(col.id))
      volumes=array(0,dim = dim(col.id))
      for (x in (1+max(c(box.size,fill.radius))):(dim(peeks)[1]-max(c(box.size,fill.radius)))){
        for (y in (1+max(c(box.size,fill.radius))):(dim(peeks)[2]-max(c(box.size,fill.radius)))){
          peeks[x,y]=(data[x,y]==max(data[(x-box.size):(x+box.size),(y-box.size):(y+box.size)]) & sum(data[(x-box.size):(x+box.size),(y-box.size):(y+box.size)]==max(data[(x-box.size):(x+box.size),(y-box.size):(y+box.size)]))==1)
          if (peeks[x,y]==TRUE){
            volumes[x,y]=vol.calc(data,x,y,fill.radius)
          }
        }
      }
      rows=c(row.id[peeks==T])
      cols=c(col.id[peeks==T])
      vols=c(volumes[peeks==T])
      spots=data.frame('x'=cols[(vols>=low.lim & vols<=high.lim)],'y'=((dim(data)[1]+1)-rows)[(vols>=low.lim & vols<=high.lim)],'row'=rows[(vols>=low.lim & vols<=high.lim)],'col'=cols[(vols>=low.lim & vols<=high.lim)],'vol'=vols[(vols>=low.lim & vols<=high.lim)])
      final=spots[order(spots$vol),]
      return(final)
    }
    if (spot.picking=='all.frames'){
      rg.calc <- function(input){
        fit=smooth.spline((1:length(input)),input,spar = 0.5)
        rang=diff(range(fit[['y']]))/sd(fit[['yin']]-fit[['y']])
        return(rang)
      }
      if (is.null(low.lim)==TRUE){
        low.lim=0
      }
      data.all=data
      data=apply(data.all,MARGIN = c(1,2),FUN = mean)
      col.id=t(array(1:dim(data)[2],dim = dim(data)[c(2,1)]))
      row.id=array(1:dim(data)[1],dim = dim(data))
      peeks=array(F,dim = dim(col.id))
      volumes=array(0,dim = dim(col.id))
      ranges=array(0,dim = dim(col.id))
      for (x in (1+max(c(box.size,fill.radius))):(dim(peeks)[1]-max(c(box.size,fill.radius)))){
        for (y in (1+max(c(box.size,fill.radius))):(dim(peeks)[2]-max(c(box.size,fill.radius)))){
          peeks[x,y]=(data[x,y]==max(data[(x-box.size):(x+box.size),(y-box.size):(y+box.size)]) & sum(data[(x-box.size):(x+box.size),(y-box.size):(y+box.size)]==max(data[(x-box.size):(x+box.size),(y-box.size):(y+box.size)]))==1)
          if (peeks[x,y]==TRUE){
            volumes[x,y]=vol.calc(data,x,y,fill.radius)
            if (volumes[x,y]>=low.lim){
              ranges[x,y]=rg.calc(apply(data.all,MARGIN = c(3),FUN = vol.calc,x=x,y=y,radius=fill.radius))
            }
          }
        }
      }
      rows=c(row.id[peeks==T])
      cols=c(col.id[peeks==T])
      vols=c(volumes[peeks==T])
      rgs=c(ranges[peeks==T])
      rg.min=3
      spots=data.frame('x'=cols[(vols>=low.lim & vols<=high.lim & rgs>=rg.min)],'y'=((dim(data)[1]+1)-rows)[(vols>=low.lim & vols<=high.lim & rgs>=rg.min)],'row'=rows[(vols>=low.lim & vols<=high.lim & rgs>=rg.min)],'col'=cols[(vols>=low.lim & vols<=high.lim & rgs>=rg.min)],'vol'=vols[(vols>=low.lim & vols<=high.lim & rgs>=rg.min)])
      final=spots[order(spots$vol),]
      return(final)
    }
  }
  #
  if (is.null(file.name)==TRUE){
    file.name=list.files(path = path.to.file,pattern = '[.]tif')
  }
  #
  image.raw=suppressWarnings(tiff::readTIFF(paste0(path.to.file,file.name),all = TRUE,as.is=TRUE))
  pixel.size=dim(image.raw[[1]])
  frame.number=length(image.raw)
  image.data=array(NA,dim = c(pixel.size[1],pixel.size[2],frame.number)); for (i in 1:frame.number){image.data[,,i]=image.raw[[i]]}
  image.avg=apply(image.data,MARGIN = c(1,2),FUN = mean)
  #
  suppressWarnings(image(t(pracma::flipud(image.avg)),col=gray.colors(cumprod(pixel.size)),x = 1:pixel.size[2],y=1:pixel.size[1],axes=FALSE,xlab='',ylab=''))
  check.1=readline(prompt = 'Image loading complete. Proceed with spot detection? (y/n):  '); if (check.1=='n'){stop('SCRIPT ABORTED BY USER')}
  #
  if (spot.picking=='composite'){
    spots=find.spots(image.avg,spot.box,spot.min,spot.max,spot.radius,spot.picking)
  }
  if (spot.picking=='all.frames'){
    spots=find.spots(image.data,spot.box,spot.min,spot.max,spot.radius,spot.picking)
  }
  points(spots$x,spots$y,col='red',cex=1.2)
  show(paste0('Median particle intensity = ',median(spots$vol)))
  show(paste0('Number of particles identified = ',length(spots$x)))
  #
  check.2=readline(prompt = 'Is initial particle selection acceptable? (y/n):  '); while (check.2=='n'){
    spot.box=as.numeric(readline(prompt = 'New spot selection box size (ENTER for default):  '))
    if (is.na(spot.box)==T){
      spot.box=6
    }
    spot.radius=as.numeric(readline(prompt = 'New spot integration radius (ENTER for default):  '))
    if (is.na(spot.radius)==T){
      spot.radius=6
    }
    spot.min=as.numeric(readline(prompt = 'New spot minimum signal (ENTER for default):  '))
    if (is.na(spot.min)==T){
      if (spot.picking=='composite'){
        spot.min=500
      }
      if (spot.picking=='all.frames'){
        spot.min=0
      }
    }
    spot.max=as.numeric(readline(prompt = 'New spot maximum signal (ENTER for default):  '))
    if (is.na(spot.max)==T){
      spot.max=Inf
    }
    if (spot.picking=='composite'){
      spots=find.spots(image.avg,spot.box,spot.min,spot.max,spot.radius,spot.picking)
    }
    if (spot.picking=='all.frames'){
      spots=find.spots(image.data,spot.box,spot.min,spot.max,spot.radius,spot.picking)
    }
    image(t(pracma::flipud(image.avg)),col=gray.colors(cumprod(pixel.size)),x = 1:pixel.size[2],y=1:pixel.size[1],axes=FALSE,xlab='',ylab='')
    points(spots$x,spots$y,col='red',cex=1.2)
    show(paste0('Median particle intensity = ',median(spots$vol)))
    show(paste0('Number of particles identified = ',length(spots$x)))
    check.2=readline(prompt = 'Is initial particle selection acceptable? (y/n/quit):  ')
    if (check.2=='quit'){
      stop('SCRIPT ABORTED BY USER')
    }
  }
  #
  particle.traces=matrix(0,nrow = frame.number,ncol = nrow(spots))
  particle.snaps=array(0,dim = c(2*spot.radius+1,2*spot.radius+1,20,nrow(spots)))
  snap.length=frame.number/20
  for (i in 1:nrow(spots)){
    particle.traces[,i]=apply(image.data[(spots$row[i]-spot.radius):(spots$row[i]+spot.radius),(spots$col[i]-spot.radius):(spots$col[i]+spot.radius),],MARGIN = c(3),FUN = sum)-apply(image.data[(spots$row[i]-spot.radius):(spots$row[i]+spot.radius),(spots$col[i]-spot.radius):(spots$col[i]+spot.radius),],MARGIN = c(3),FUN = median)*((2*spot.radius+1)^2)
    for (j in 1:20){
      particle.snaps[,,j,i]=apply(image.data[(spots$row[i]-spot.radius):(spots$row[i]+spot.radius),(spots$col[i]-spot.radius):(spots$col[i]+spot.radius),((j-1)*snap.length+1):(j*snap.length)],MARGIN = c(1,2),FUN = mean)
    }
  }
  row.names(particle.traces)=time.step*(1:frame.number)
  #
  dev.print(pdf,paste0(path.to.file,'Identified_Spots.pdf'))
  #
  utils::write.table(spots,file = paste0(path.to.file,'initial_particle_summary.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  utils::write.table(particle.traces,file = paste0(path.to.file,'initial_particle_traces.txt'),quote = FALSE,sep = '\t',row.names = TRUE,col.names = FALSE)
  save(list = c('image.avg','spots','particle.traces','particle.snaps','time.step','frame.number','pixel.size','spot.radius'),file = paste0(path.to.file,'Initial-Particle_Data.RData'))
  id.spots.output.files=list('composite image data'=image.avg,'initial particle summary'=spots,'particle traces'=particle.traces,'particle images data'=particle.snaps)
  return(id.spots.output.files)
}

refine.particles <- function(path.to.file='./',file.name='Initial-Particle_Data.RData',skip.manual='n',signal.step=NULL,auto.filter='none',classification.strategy='classic',background.subtraction='lower.quartile'){
  count.condense <- function(input){
    new.count=input
    for (i in 1:length(input)){
      new.count[i]=length(unique(input[1:i]))
    }
    return(new.count)
  }
  lq.calc <- function(input){
    input.2=na.omit(c(input))
    result=(input.2[order(input.2)])[round(0.25*length(input.2))]
    return(result)
  }
  k.mode <- function(data){
    k=density(data)
    mode=k[['x']][which.max(k[['y']])]
    return(mode)
  }
  classify.states <- function(data,signal.step,fun='classic',background=0){
    if (fun=='classic'){
      fit=smooth.spline((1:length(data))[is.na(data)==F],na.omit(data),spar = 0.5)
      d1=predict(object = fit,x = 1:length(data),deriv = 1)
      swaps=c(1,d1[['x']][(abs(d1[['y']])>=1.0) & (splus2R::peaks(x = abs(d1[['y']]),span = 5))],max(d1[['x']]))
      pos.averages=rep(NA,times=length(data))
      states=rep(0,times=length(data))
      state.averages=rep(NA,times=length(data))
      for (i in 1:(length(swaps)-1)){
        pos.averages[swaps[i]:swaps[i+1]]=mean(na.omit(data[swaps[i]:swaps[i+1]]))
      }
      for (j in 0:50){
        states[(pos.averages>=(signal.step*j-signal.step/2+background)) & (pos.averages<=(signal.step*j+signal.step/2+background))]=j
      }
      for (k in 0:50){
        state.averages[states==k]=mean(na.omit(data[states==k]))
      }
      final=as.data.frame(cbind(states,state.averages))
      colnames(final)=c('id','avg')
      return(final)
    }
    if (fun=='basic'){
      pos.averages=matrix(NA,ncol = ncol(data),nrow = nrow(data))
      for (n in 1:ncol(data)){
        fit=smooth.spline((1:length(data[,n]))[is.na(data[,n])==F],na.omit(data[,n]),spar = 0.5)
        d1=predict(object = fit,x = 1:length(data[,n]),deriv = 1)
        swaps=c(1,d1[['x']][(abs(d1[['y']])>=1.0) & (splus2R::peaks(x = abs(d1[['y']]),span = 5))],max(d1[['x']]))
        for (i in 1:(length(swaps)-1)){
          pos.averages[swaps[i]:swaps[i+1],n]=mean(na.omit(data[swaps[i]:swaps[i+1],n]))
        }
      }
      final=pos.averages
      return(final)
    }
    if (fun=='uni.dbscan'){
      pos.averages=matrix(NA,ncol = ncol(data),nrow = nrow(data))
      for (n in 1:ncol(data)){
        fit=smooth.spline((1:length(data[,n]))[is.na(data[,n])==F],na.omit(data[,n]),spar = 0.5)
        d1=predict(object = fit,x = 1:length(data[,n]),deriv = 1)
        swaps=c(1,d1[['x']][(abs(d1[['y']])>=1.0) & (splus2R::peaks(x = abs(d1[['y']]),span = 5))],max(d1[['x']]))
        for (i in 1:(length(swaps)-1)){
          pos.averages[swaps[i]:swaps[i+1],n]=mean(na.omit(data[swaps[i]:swaps[i+1],n]))
        }
      }
      uni.dbscan <- function(datter,eps){
        temp.0=order(c(datter))
        temp.1=datter[temp.0]
        temp.2=c(abs(temp.1[1:(length(temp.1)-1)]-temp.1[2:length(temp.1)]),0)>eps
        temp.3=c(1,(1:length(temp.2))[temp.2],length(temp.2))
        temp.4=rep(NA,times=length(temp.2))
        for (i in 0:(length(temp.3)-2)){
          temp.4[(temp.3[i+1]):(temp.3[i+2])]=i
        }
        temp.5=matrix(temp.4[order(temp.0)],nrow = nrow(datter),ncol = ncol(datter))
        return(temp.5)
      }
      classes=uni.dbscan(pos.averages,eps = 1*sd(na.omit(data-pos.averages)))
      values=matrix(NA,nrow = nrow(classes),ncol = ncol(classes))
      for (k in 0:max(classes)){
        values[classes==k]=mean(na.omit(data[classes==k]))
      }
      final=list('id'=classes,'avg'=values)
      return(final)
    }
    if (fun=='var.shift'){
      pos.averages=rep(NA,times=length(data))
      fit=smooth.spline((1:length(data))[is.na(data)==F],na.omit(data),spar = 0.5)
      d1=predict(object = fit,x = 1:length(data),deriv = 1)
      swaps=c(1,d1[['x']][(abs(d1[['y']])>=1.0) & (splus2R::peaks(x = abs(d1[['y']]),span = 5))],max(d1[['x']]))
      for (i in 1:(length(swaps)-1)){
        pos.averages[swaps[i]:swaps[i+1]]=mean(na.omit(data[swaps[i]:swaps[i+1]]))
      }
      var.shift <- function(datter,eps){
        temp.0=order(datter)
        temp.1=datter[temp.0]
        temp.2=c(abs(temp.1[1:(length(temp.1)-1)]-temp.1[2:length(temp.1)]),0)>eps
        temp.3=c(1,(1:length(temp.2))[temp.2],length(temp.2))
        temp.4=rep(NA,times=length(temp.2))
        for (i in 1:(length(temp.3)-1)){
          temp.4[(temp.3[i]):(temp.3[i+1])]=i
        }
        temp.5=temp.4[order(temp.0)]
        return(temp.5)
      }
      classes=var.shift(pos.averages,eps = 1.7*sd(na.omit(data-pos.averages)))
      values=rep(NA,times=length(classes))
      for (k in 0:max(classes)){
        values[classes==k]=mean(na.omit(data[classes==k]))
      }
      classes[abs(values)<=2*sd(na.omit(data-pos.averages))]=0
      final=data.frame('id'=classes,'avg'=values)
      return(final)
    }
  }
  load(file = paste0(path.to.file,file.name))
  #
  particle.trace.rolls=apply(particle.traces,MARGIN = c(2),FUN = data.table::frollmean,n=5,align='center')
  if (background.subtraction=='k.mode'){
    background=k.mode(na.omit(particle.trace.rolls))
    particle.trace.rolls=particle.trace.rolls-background
  }
  if (background.subtraction=='basal.states'){
    background=median(na.omit(apply(classify.states(particle.trace.rolls,fun = 'basic'),MARGIN = c(2),FUN = min)))
    particle.trace.rolls=particle.trace.rolls-background
  }
  if (background.subtraction=='lower.quartile'){
    background=lq.calc(classify.states(particle.trace.rolls,fun = 'basic'))
    particle.trace.rolls=particle.trace.rolls-background
  }
  particles.to.keep=rep(FALSE,times=nrow(spots))
  residence.times=c('Particle','Start','Stop','Residence','State Signal')
  dwell.calls=c('Particle','State','Dwell')
  residence.calls=c('Particle','Bound','Residence')
  if (classification.strategy=='classic' | classification.strategy=='var.shift'){
    state.calls=matrix(NA,nrow = frame.number,ncol = nrow(spots))
    if (is.null(signal.step)==TRUE){
      signal.step=3*sd(na.omit(particle.trace.rolls-classify.states(particle.trace.rolls,fun = 'basic')))
    }
  }
  if (classification.strategy=='uni.dbscan'){
    classifications=classify.states(particle.trace.rolls,fun = 'uni.dbscan')
    state.calls=classifications[['id']]
  }
  filtering=rep(F,times=nrow(spots))
  COUNTER=0
  for (i in 1:nrow(spots)){
    if (classification.strategy=='classic' | classification.strategy=='var.shift'){
      states=classify.states(particle.trace.rolls[,i],signal.step,fun = classification.strategy)
      state.calls[,i]=states$id
    }
    dwell.calls=rbind(dwell.calls,cbind(rep(i,times=length(rle(state.calls[,i])[['values']])),rle(state.calls[,i])[['values']],rle(state.calls[,i])[['lengths']]*time.step))
    residence.calls=rbind(residence.calls,cbind(rep(i,times=length(rle(state.calls[,i]>=1)[['values']])),rle(state.calls[,i]>=1)[['values']],rle(state.calls[,i]>=1)[['lengths']]*time.step))
    if  (auto.filter=='unbound') {
      filtering[i]=(sum(state.calls[,i])==0)
    }
    if  (auto.filter=='all.stable') {
      filtering[i]=(length(table(state.calls[,i]))==1)
    }
  }
  #
  plot((1:frame.number)*time.step,apply(state.calls>0,MARGIN = c(1),FUN = sum)/nrow(spots),type='l',main='Photobleaching Check',xlab='Time (s)',ylab='Proportion of Bound Particles')
  temp.3=readline('Proceed with particle refinement? (y/n):  ')
  if (temp.3=='n'){
    stop('SCRIPT ABORTED BY USER')
  }
  #
  for (i in 1:nrow(spots)){
    if (skip.manual=='n' & filtering[i]==FALSE){
      if (classification.strategy=='classic' | classification.strategy=='var.shift'){
        states=classify.states(particle.trace.rolls[,i],signal.step,fun = classification.strategy)
      }
      if (classification.strategy=='uni.dbscan'){
        states=data.frame('id'=classifications[['id']][,i],'avg'=classifications[['avg']][,i])
      }
      par(fig=c(0,1,0.6,1))
      temp.1=matrix(particle.snaps[,,,i],nrow = (2*spot.radius+1),ncol = (2*spot.radius+1)*20)
      image(t(pracma::flipud(rbind(temp.1[,1:(ncol(temp.1)/2)],temp.1[,(ncol(temp.1)/2+1):ncol(temp.1)]))),col=gray.colors(length(temp.1)),axes=FALSE)
      par(fig=c(0,1,0,0.75),new = TRUE)
      plot((1:frame.number)*time.step,particle.trace.rolls[,i],type='l',main = 'Particle Trace',xlab='Time (s)',ylab = 'Signal')
      lines(((1:frame.number)*time.step)[is.na(states$avg)==F],states$avg[is.na(states$avg)==F],col='green',lwd=3)
      temp.2=readline('Should this particle be used for analysis? (y/n/quit/finish):  ')
      if (temp.2=='y'){
        particles.to.keep[i]=TRUE
        COUNTER=COUNTER+1
        event.number=as.numeric(readline('How many binding events will you record for this particle?:  '))
        if (length(event.number)==1 & is.na(event.number)==FALSE & event.number>0){
          event.times=matrix(0,nrow=event.number,ncol = 5); event.times[,1]=rep(COUNTER,times=event.number)
          for (j in 1:event.number){
            show('Please click on the starting point then stopping point of a binding event (bottom graph -- green fit line).')
            event.times[j,2:3]=as.numeric(identify(((1:frame.number)*time.step)[is.na(states$avg)==F],states$avg[is.na(states$avg)==F],n=2)*time.step)
            event.times[j,4]=diff(event.times[j,2:3])
            event.times[j,5]=mean(na.omit(particle.trace.rolls[(event.times[j,2]/time.step):(event.times[j,3]/time.step),i]))
            arrows(x0 = event.times[j,2],x1 = event.times[j,3],y0=mean(na.omit(particle.trace.rolls[(event.times[j,2]/time.step):(event.times[j,3]/time.step),i])),y1=mean(na.omit(particle.trace.rolls[(event.times[j,2]/time.step):(event.times[j,3]/time.step),i])),angle = 90,code = 3,col = 'red')
          }
          residence.times=rbind(residence.times,event.times)
        }
        dev.print(pdf,paste0(path.to.file,'Particle_Trace_',COUNTER,'.pdf'))
      }
      if (temp.2=='quit'){
        stop('SCRIPT ABORTED BY USER')
      }
    }
    if (skip.manual=='y' & i==1){
      temp.2='blank'
    }
    if (temp.2=='finish'){
      break
    }
  }
  show(paste0('Particle refinement completed -- beginning data exports!'))
  residence.calls=as.data.frame(as.matrix(residence.calls[2:nrow(residence.calls),])); for (k in 1:3){residence.calls[,k]=as.numeric(residence.calls[,k])}; colnames(residence.calls)=c('Particle','Bound','Residence')
  dwell.calls=as.data.frame(as.matrix(dwell.calls[2:nrow(dwell.calls),])); for (k in 1:3){dwell.calls[,k]=as.numeric(dwell.calls[,k])}; colnames(dwell.calls)=c('Particle','State','Dwell')
  if (skip.manual=='n'){
    if (length(dim(residence.times))==2){
      residence.data=as.data.frame(as.matrix(residence.times[2:nrow(residence.times),])); for (k in 1:4){residence.data[,k]=as.numeric(residence.data[,k])}; colnames(residence.data)=c('Particle','Start','Stop','Residence')
    } else {
      residence.data=c('Particle','Start','Stop','Residence')
    }
    refined.particle.traces=particle.traces[,particles.to.keep]
    refined.particle.trace.rolls=particle.trace.rolls[,particles.to.keep]
    refined.spots=spots[particles.to.keep,]
    refined.particle.snaps=particle.snaps[,,,particles.to.keep]
    refined.state.calls=state.calls[,particles.to.keep]
    refined.residence.calls=residence.calls[which(residence.calls[,1]%in%which(particles.to.keep)),]; refined.residence.calls[,1]=count.condense(refined.residence.calls[,1])
    refined.dwell.calls=dwell.calls[which(dwell.calls[,1]%in%which(particles.to.keep)),]; refined.dwell.calls[,1]=count.condense(refined.dwell.calls[,1])
  }
  #
  if (skip.manual=='n'){
    save(list = c('image.avg','residence.data','refined.particle.traces','refined.particle.trace.rolls','refined.spots','refined.particle.snaps','state.calls','residence.calls','dwell.calls','refined.state.calls','refined.residence.calls','refined.dwell.calls'),file = paste0(path.to.file,'Refined-Particle_Data.RData'))
    utils::write.table(residence.data,file = paste0(path.to.file,'residence_data.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
    utils::write.table(refined.particle.traces,file = paste0(path.to.file,'selected_particle_traces.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
    utils::write.table(refined.particle.trace.rolls,file = paste0(path.to.file,'selected_particle_smoothed-traces.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
    utils::write.table(refined.spots,file = paste0(path.to.file,'selected_particle_summary.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
    utils::write.table(refined.state.calls,file = paste0(path.to.file,'selected_particle_state-calls.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
    utils::write.table(refined.residence.calls,file = paste0(path.to.file,'selected_particle_residence-calls.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
    utils::write.table(refined.dwell.calls,file = paste0(path.to.file,'selected_particle_dwell-calls.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
    return(list('image.avg'=image.avg,'residence.data'=residence.data,'refined.particle.traces'=refined.particle.traces,'refined.particle.trace.rolls'=refined.particle.trace.rolls,'refined.spots'=refined.spots,'refined.particle.snaps'=refined.particle.snaps,'state.calls'=state.calls,'residence.calls'=residence.calls,'dwell.calls'=dwell.calls))
  }
  if (skip.manual=='y'){
    save(list = c('image.avg','state.calls','residence.calls','dwell.calls'),file = paste0(path.to.file,'Refined-Particle_Data.RData'))
    return(list('image.avg'=image.avg,'state.calls'=state.calls,'residence.calls'=residence.calls,'dwell.calls'=dwell.calls))
  }
  utils::write.table(state.calls,file = paste0(path.to.file,'all-particle_state-calls.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  utils::write.table(dwell.calls,file = paste0(path.to.file,'all-particle_dwell-calls.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  utils::write.table(residence.calls,file = paste0(path.to.file,'all-particle_residence-calls.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
}

calc.kn1 <- function(path.to.file='./',file.name='Refined-Particle_Data.RData',use.auto.times='n',min.residence=0,max.residence=Inf){
  load(file = paste0(path.to.file,file.name))
  #
  if (use.auto.times=='n'){
    residence.times=residence.data$Residence
  }
  if (use.auto.times=='y'){
    residence.times=na.omit(refined.residence.calls$Residence[refined.residence.calls$Bound==T])
  }
  #
  residence.time.values=residence.times[((residence.times>=min.residence) & (residence.times<=max.residence))]
  residence.time.outliers=residence.times[((residence.times<min.residence) | (residence.times>max.residence))]
  #
  mod1=fitdistrplus::fitdist(data=as.numeric(residence.time.values),distr = 'gamma',method = c('mle'))
  mod1.shape=mod1$estimate[1]
  mod1.rate=mod1$estimate[2]
  #
  fit.x=seq(0,1.5*max(residence.time.values),0.01)
  mod1.fit=dgamma(fit.x,shape = mod1.shape,rate = mod1.rate)
  #
  par(fig =c(0,1,0.3,1))
  plot(fit.x,mod1.fit,type='l',xlab = 'Residence Time (s)',ylab = 'Probability Density',main = 'Distribution of Residence Times',col='blue')
  lines(density(residence.time.values),col='red')
  legend('topright',legend = c('Kernel Density','Gamma Fit','Data Points','Outliers'),fill = c('red','blue','black','grey'),col= c('red','blue','black','grey'))
  text(0.5*median(fit.x),max(mod1.fit),labels = paste0('Average Residence Time = ',signif(mean(residence.time.values),3),' (s)'),adj=c(0,1))
  text(0.5*median(fit.x),max(mod1.fit),labels = paste0('Dissociation Rate = ',signif(1/mean(residence.time.values),3),' (1/s)'),adj=c(0,3))
  #
  par(fig =c(0,1,0,0.4),new = TRUE)
  plot(residence.time.values,rep(0,times=length(residence.time.values)),type = 'p',col='black',main = NULL,xlim = range(fit.x),ann=FALSE,yaxt='n',xaxt='n')
  points(residence.time.outliers,rep(0,times=length(residence.time.outliers)),col='grey')
  abline(v=mean(residence.time.values),col='green',lwd=3)
  #
  show(paste0('Average Residence Time = ',signif(mean(residence.time.values),3),' s'))
  show(paste0('Dissociation Rate = ',signif(1/mean(residence.time.values),3),' 1/s'))
  #
  check.1=readline('Do you wish to save this analysis? (y/n):   ')
  if (check.1=='y'){
    dev.print(pdf,paste0(path.to.file,'Dissociation-Rate_Analysis.pdf'))
  }
}

FRET.align <- function(path.to.file='./',file.name=NULL,alignment=list('r'=0.97,'theta'=7),search.radius=6,integration.radius=6,spot.min=NULL,spot.max=Inf){
  cart2pol <- function(data,centroid){
    data.pol=data.frame('r'=sqrt((data$x-centroid[1])^2+(data$y-centroid[2])^2),'theta'=atan2((data$y-centroid[2]),(data$x-centroid[1]))*180/pi)
    return(data.pol)
  }
  pol2cart <- function(data,centroid){
    data.cart=data.frame('x'=data$r*cos(data$theta*pi/180)+centroid[1],'y'=data$r*sin(data$theta*pi/180)+centroid[2])
    return(data.cart)
  }
  cart2mat <- function(data,row.num){
    data.mat=data.frame('row'=round(row.num-data$y+1),'col'=round(data$x))
    return(data.mat)
  }
  mat2cart <- function(data,row.num){
    data.cart=data.frame('x'=data$col,'y'=row.num-data$row+1)
    return(data.cart)
  }
  pol.adjust <- function(data,r.factor,theta.factor){
    data.adj=data.frame('r'=data$r*r.factor,'theta'=data$theta+theta.factor)
    return(data.adj)
  }
  get.msd <- function(par,Cy5.spots,Cy3.scale){
    cart2pol <- function(data,centroid){
      data.pol=data.frame('r'=sqrt((data$x-centroid[1])^2+(data$y-centroid[2])^2),'theta'=atan2((data$y-centroid[2]),(data$x-centroid[1]))*180/pi)
      return(data.pol)
    }
    pol2cart <- function(data,centroid){
      data.cart=data.frame('x'=data$r*cos(data$theta*pi/180)+centroid[1],'y'=data$r*sin(data$theta*pi/180)+centroid[2])
      return(data.cart)
    }
    cart2mat <- function(data,row.num){
      data.mat=data.frame('row'=round(row.num-data$y+1),'col'=round(data$x))
      return(data.mat)
    }
    mat2cart <- function(data,row.num){
      data.cart=data.frame('x'=data$col,'y'=row.num-data$row+1)
      return(data.cart)
    }
    pol.adjust <- function(data,r.factor,theta.factor){
      data.adj=data.frame('r'=data$r*r.factor,'theta'=data$theta+theta.factor)
      return(data.adj)
    }
    mat.id.grab <- function(row.col,data){
      value=data[row.col[1],row.col[2]]
      return(value)
    }
    row.num=dim(Cy3.scale)[1]
    centroid=rev(dim(Cy3.scale)/2+0.5)
    r.factor=par[1]
    theta.factor=par[2]
    #
    Cy3.projections=cart2mat(pol2cart(pol.adjust(cart2pol(Cy5.spots[,1:2],centroid),r.factor,theta.factor),centroid),row.num)
    Cy3.projections.filtered=Cy3.projections[(Cy3.projections$row<=dim(Cy3.scale)[1] & Cy3.projections$col<=dim(Cy3.scale)[2]),]
    Cy3.projection.signals=apply(X = Cy3.projections.filtered,MARGIN = 1,FUN = mat.id.grab,data=Cy3.scale)
    Cy5.signals.filtered=Cy5.spots[(Cy3.projections$row<=dim(Cy3.scale)[1] & Cy3.projections$col<=dim(Cy3.scale)[2]),5]
    msd=mean((Cy3.projection.signals-Cy5.signals.filtered)^2)
    return(msd)
  }
  woh.scale <- function(data){
    a=c(data)
    b=(data-min(a))/diff(range(a))
    return(b)
  }
  auto.align <- function(Cy5.spots,Cy3.scale,pars){
    get.msd <- function(par,Cy5.spots,Cy3.scale){
      cart2pol <- function(data,centroid){
        data.pol=data.frame('r'=sqrt((data$x-centroid[1])^2+(data$y-centroid[2])^2),'theta'=atan2((data$y-centroid[2]),(data$x-centroid[1]))*180/pi)
        return(data.pol)
      }
      pol2cart <- function(data,centroid){
        data.cart=data.frame('x'=data$r*cos(data$theta*pi/180)+centroid[1],'y'=data$r*sin(data$theta*pi/180)+centroid[2])
        return(data.cart)
      }
      cart2mat <- function(data,row.num){
        data.mat=data.frame('row'=round(row.num-data$y+1),'col'=round(data$x))
        return(data.mat)
      }
      mat2cart <- function(data,row.num){
        data.cart=data.frame('x'=data$col,'y'=row.num-data$row+1)
        return(data.cart)
      }
      pol.adjust <- function(data,r.factor,theta.factor){
        data.adj=data.frame('r'=data$r*r.factor,'theta'=data$theta+theta.factor)
        return(data.adj)
      }
      mat.id.grab <- function(row.col,data){
        value=data[row.col[1],row.col[2]]
        return(value)
      }
      row.num=dim(Cy3.scale)[1]
      centroid=rev(dim(Cy3.scale)/2+0.5)
      r.factor=par[1]
      theta.factor=par[2]
      #
      Cy3.projections=cart2mat(pol2cart(pol.adjust(cart2pol(Cy5.spots[,1:2],centroid),r.factor,theta.factor),centroid),row.num)
      Cy3.projections.filtered=Cy3.projections[(Cy3.projections$row<=dim(Cy3.scale)[1] & Cy3.projections$col<=dim(Cy3.scale)[2]),]
      Cy3.projection.signals=apply(X = Cy3.projections.filtered,MARGIN = 1,FUN = mat.id.grab,data=Cy3.scale)
      Cy5.signals.filtered=Cy5.spots[(Cy3.projections$row<=dim(Cy3.scale)[1] & Cy3.projections$col<=dim(Cy3.scale)[2]),5]
      msd=mean((Cy3.projection.signals-Cy5.signals.filtered)^2)
      return(msd)
    }
    model.fit=optim(par = c(pars[['r']],pars[['theta']]),fn = get.msd,method = "L-BFGS-B",lower = c(0.8,-45),upper = c(1.09,45),Cy5.spots = Cy5.spots,Cy3.scale = Cy3.scale,control = list('ndeps'=c(0.01,0.5)))
    fit.parameters=list('r'=model.fit[['par']][1],'theta'=model.fit[['par']][2])
    return(fit.parameters)
  }
  lq.calc <- function(input){
    input.2=na.omit(c(input))
    result=(input.2[order(input.2)])[round(0.25*length(input.2))]
    return(result)
  }
  find.spots <- function(data,box.size,low.lim,high.lim,fill.radius){
    lq.calc <- function(input){
      input.2=na.omit(c(input))
      result=(input.2[order(input.2)])[round(0.25*length(input.2))]
      return(result)
    }
    vol.calc <- function(input,x,y,radius){
      lq.calc <- function(input){
        input.2=na.omit(c(input))
        result=(input.2[order(input.2)])[round(0.25*length(input.2))]
        return(result)
      }
      result=sum(input[(x-radius):(x+radius),(y-radius):(y+radius)])-lq.calc(input[(x-radius):(x+radius),(y-radius):(y+radius)])*(2*radius+1)^2
      return(result)
    }
    col.id=t(array(1:dim(data)[2],dim = dim(data)[c(2,1)]))
    row.id=array(1:dim(data)[1],dim = dim(data))
    peeks=array(F,dim = dim(col.id))
    volumes=array(0,dim = dim(col.id))
    for (x in (1+max(c(box.size,fill.radius))):(dim(peeks)[1]-max(c(box.size,fill.radius)))){
      for (y in (1+max(c(box.size,fill.radius))):(dim(peeks)[2]-max(c(box.size,fill.radius)))){
        peeks[x,y]=(data[x,y]==max(data[(x-box.size):(x+box.size),(y-box.size):(y+box.size)]) & sum(data[(x-box.size):(x+box.size),(y-box.size):(y+box.size)]==max(data[(x-box.size):(x+box.size),(y-box.size):(y+box.size)]))==1)
        if (peeks[x,y]==TRUE){
          volumes[x,y]=vol.calc(data,x,y,fill.radius)
        }
      }
    }
    rows=c(row.id[peeks==T])
    cols=c(col.id[peeks==T])
    vols=c(volumes[peeks==T])
    spots=data.frame('x'=cols[(vols>=low.lim & vols<=high.lim)],'y'=((dim(data)[1]+1)-rows)[(vols>=low.lim & vols<=high.lim)],'row'=rows[(vols>=low.lim & vols<=high.lim)],'col'=cols[(vols>=low.lim & vols<=high.lim)],'vol'=vols[(vols>=low.lim & vols<=high.lim)])
    final=spots[order(spots$vol),]
    return(final)
  }
  #
  if (is.null(file.name)==TRUE){
    file.name=list.files(path = path.to.file,pattern = '[.]tif')
    show(paste0('Cy3 Filename = ',file.name[1]))
    show(paste0('Cy5 Filename = ',file.name[2]))
  }
  image.raw.Cy3=suppressWarnings(tiff::readTIFF(paste0(path.to.file,file.name[1]),all = TRUE,as.is=TRUE))
  frame.number.Cy3=length(image.raw.Cy3)
  pixel.size.Cy3=dim(image.raw.Cy3[[1]])
  image.data.Cy3=array(NA,dim = c(pixel.size.Cy3[1],pixel.size.Cy3[2],frame.number.Cy3)); for (i in 1:frame.number.Cy3){image.data.Cy3[,,i]=image.raw.Cy3[[i]]}
  Cy3.ref=apply(image.data.Cy3,MARGIN = c(1,2),FUN = mean)
  Cy3.scale=woh.scale(Cy3.ref)
  image.raw.Cy5=suppressWarnings(tiff::readTIFF(paste0(path.to.file,file.name[2]),all = TRUE,as.is=TRUE))
  frame.number.Cy5=length(image.raw.Cy5)
  pixel.size.Cy5=dim(image.raw.Cy5[[1]])
  image.data.Cy5=array(NA,dim = c(pixel.size.Cy5[1],pixel.size.Cy5[2],frame.number.Cy5)); for (i in 1:frame.number.Cy5){image.data.Cy5[,,i]=image.raw.Cy5[[i]]}
  Cy5.ref=apply(image.data.Cy5,MARGIN = c(1,2),FUN = mean)
  Cy5.scale=woh.scale(Cy5.ref)
  par(fig = c(0,0.5,0,1))
  suppressWarnings(image(t(pracma::flipud(Cy3.ref)),col=gray.colors(cumprod(pixel.size.Cy3)),x = 1:pixel.size.Cy3[2],y=1:pixel.size.Cy3[1],axes=FALSE,xlab='',ylab='',main = 'Cy3'))
  par(fig = c(0.5,1,0,1),new  = TRUE)
  suppressWarnings(image(t(pracma::flipud(Cy5.ref)),col=gray.colors(cumprod(pixel.size.Cy5)),x = 1:pixel.size.Cy5[2],y=1:pixel.size.Cy5[1],axes=FALSE,xlab='',ylab='',main = 'Cy5'))
  check.1=readline(prompt = 'Image loading complete. Proceed with alignment? (y/n):  '); if (check.1=='n'){stop('SCRIPT ABORTED BY USER')}
  #
  if (is.null(spot.min)==TRUE){
    spot.min=5*lq.calc(Cy3.scale)*(2*integration.radius+1)^2
  }
  Cy5.spots=find.spots(Cy5.scale,search.radius,spot.min,spot.max,integration.radius)
  par(fig = c(0.25,0.75,0,1))
  suppressWarnings(image(t(pracma::flipud(Cy5.ref)),col=gray.colors(cumprod(pixel.size.Cy5)),x = 1:pixel.size.Cy5[2],y=1:pixel.size.Cy5[1],axes=FALSE,xlab='',ylab='',main = 'Cy5 Reference Points'))
  points(Cy5.spots$x,Cy5.spots$y,col = 'red')
  show(paste0('Median spot volume = ',median(Cy5.spots$vol)))
  show(paste0('Number of spots identified = ',length(Cy5.spots$x)))
  check.2=readline(prompt = 'Is reference spot selection acceptable? (y/n):  '); if (check.2!='n' & check.2!='y'){stop('INVALID RESPONSE')}
  while (check.2=='n'){
    search.radius=as.numeric(readline(prompt = 'New search box radius (ENTER for default):  '))
    if (is.na(search.radius)==T){
      search.radius=6
    }
    integration.radius=as.numeric(readline(prompt = 'New integration box radius (ENTER for default):  '))
    if (is.na(integration.radius)==T){
      integration.radius=6
    }
    spot.min=as.numeric(readline(prompt = 'New spot minimum volume (ENTER for default):  '))
    if (is.na(spot.min)==T){
      spot.min=5*lq.calc(Cy3.scale)*(2*integration.radius+1)^2
    }
    spot.max=as.numeric(readline(prompt = 'New spot maximum volume (ENTER for default):  '))
    if (is.na(spot.max)==T){
      spot.max=Inf
    }
    Cy5.spots=find.spots(Cy5.scale,search.radius,spot.min,spot.max,integration.radius)
    par(fig = c(0.25,0.75,0,1))
    suppressWarnings(image(t(pracma::flipud(Cy5.ref)),col=gray.colors(cumprod(pixel.size.Cy5)),x = 1:pixel.size.Cy5[2],y=1:pixel.size.Cy5[1],axes=FALSE,xlab='',ylab='',main = 'Cy5 Reference Points'))
    points(Cy5.spots$x,Cy5.spots$y,col = 'red')
    show(paste0('Median spot volume = ',median(Cy5.spots$vol)))
    show(paste0('Number of spots identified = ',length(Cy5.spots$x)))
    check.2=readline(prompt = 'Is reference spot selection acceptable? (y/n/quit):  '); if (check.2!='n' & check.2!='y' & check.2!='quit'){stop('INVALID RESPONSE')}; if (check.2=='quit'){stop('SCRIPT ABORTED BY USER')}
  }
  Cy3.spots.matched=find.spots(Cy3.scale,search.radius,spot.min,spot.max,integration.radius)
  #
  if (is.list(alignment)==TRUE){
    auto.pars=auto.align(Cy5.spots,Cy3.scale,alignment)
  }
  if (is.list(alignment)==FALSE & alignment=='manual'){
    auto.pars=list('r'=1,'theta'=0)
  }
  r.adj.value=auto.pars[['r']]
  theta.adj.value=auto.pars[['theta']]
  #
  par(fig = c(0,0.5,0,1))
  suppressWarnings(image(t(pracma::flipud(Cy3.ref)),col=gray.colors(cumprod(pixel.size.Cy3)),x = 1:pixel.size.Cy3[2],y=1:pixel.size.Cy3[1],axes=FALSE,xlab='',ylab='',main = 'Cy3'))
  points(Cy5.spots$x,Cy5.spots$y,col = 'red',pch = 'x',cex=0.4)
  points(pol2cart(pol.adjust(cart2pol(Cy5.spots[,1:2],rev(dim(Cy3.scale)/2+0.5)),r.adj.value,theta.adj.value),rev(dim(Cy3.scale)/2+0.5)),col = 'green')
  par(fig = c(0.5,1,0,1),new  = TRUE)
  suppressWarnings(image(t(pracma::flipud(Cy5.ref)),col=gray.colors(cumprod(pixel.size.Cy5)),x = 1:pixel.size.Cy5[2],y=1:pixel.size.Cy5[1],axes=FALSE,xlab='',ylab='',main = 'Cy5'))
  points(pol2cart(pol.adjust(cart2pol(Cy5.spots[,1:2],rev(dim(Cy3.scale)/2+0.5)),r.adj.value,theta.adj.value),rev(dim(Cy3.scale)/2+0.5)),col = 'green',pch = 'x',cex=0.4)
  points(Cy5.spots$x,Cy5.spots$y,col = 'red')
  show(paste0('r = ',signif(r.adj.value,3)))
  show(paste0('theta = ',signif(theta.adj.value,3)))
  show(paste0('RMSD = ',signif(sqrt(get.msd(c(r.adj.value,theta.adj.value),Cy5.spots,Cy3.scale)),3)))
  check.3=readline(prompt = 'Is initial alignment acceptable? (y/n):  '); if (check.3!='n' & check.3!='y'){stop('INVALID RESPONSE')}
  while (check.3=='n'){
    r.adj.value=as.numeric(readline(prompt = 'New scaling factor (ENTER for default):  '))
    if (is.na(r.adj.value)==T){
      r.adj.value=auto.pars[['r']]
    }
    theta.adj.value=as.numeric(readline(prompt = 'New rotational adjustment factor (ENTER for default):  '))
    if (is.na(theta.adj.value)==T){
      theta.adj.value=auto.pars[['theta']]
    }
    par(fig = c(0,0.5,0,1))
    suppressWarnings(image(t(pracma::flipud(Cy3.ref)),col=gray.colors(cumprod(pixel.size.Cy3)),x = 1:pixel.size.Cy3[2],y=1:pixel.size.Cy3[1],axes=FALSE,xlab='',ylab='',main = 'Cy3'))
    points(Cy5.spots$x,Cy5.spots$y,col = 'red',pch = 'x',cex=0.4)
    points(pol2cart(pol.adjust(cart2pol(Cy5.spots[,1:2],rev(dim(Cy3.scale)/2+0.5)),r.adj.value,theta.adj.value),rev(dim(Cy3.scale)/2+0.5)),col = 'green')
    par(fig = c(0.5,1,0,1),new  = TRUE)
    suppressWarnings(image(t(pracma::flipud(Cy5.ref)),col=gray.colors(cumprod(pixel.size.Cy5)),x = 1:pixel.size.Cy5[2],y=1:pixel.size.Cy5[1],axes=FALSE,xlab='',ylab='',main = 'Cy5'))
    points(pol2cart(pol.adjust(cart2pol(Cy5.spots[,1:2],rev(dim(Cy3.scale)/2+0.5)),r.adj.value,theta.adj.value),rev(dim(Cy3.scale)/2+0.5)),col = 'green',pch = 'x',cex=0.4)
    points(Cy5.spots$x,Cy5.spots$y,col = 'red')
    show(paste0('r = ',signif(r.adj.value,3)))
    show(paste0('theta = ',signif(theta.adj.value,3)))
    show(paste0('RMSD = ',signif(sqrt(get.msd(c(r.adj.value,theta.adj.value),Cy5.spots,Cy3.scale)),3)))
    check.3=readline(prompt = 'Is alignment acceptable? (y/n/quit):  '); if (check.3=='quit'){stop('SCRIPT ABORTED BY USER')}; if (check.3!='n' & check.3!='y' & check.3!='quit'){stop('INVALID RESPONSE')}
  }
  # 
  save(r.adj.value,theta.adj.value,file = paste0(path.to.file,'Alignment_Parameters.RData'))
}

FRET.id <- function(time.step,path.to.file='./',Cy3.file.name=NULL,Cy5.file.name=NULL,align.file.name='Alignment_Parameters.RData',Cy3.search.box=6,Cy5.search.box=6,Cy3.integration.radius=6,Cy5.integration.radius=6,Cy3.spot.min=3500,Cy5.spot.min=1e4,Cy3.spot.max=Inf,Cy5.spot.max=Inf){
  lq.calc <- function(input){
    input.2=na.omit(c(input))
    result=(input.2[order(input.2)])[round(0.25*length(input.2))]
    return(result)
  }
  vol.calc <- function(input,row,col,radius){
    lq.calc <- function(input){
      input.2=na.omit(c(input))
      result=(input.2[order(input.2)])[round(0.25*length(input.2))]
      return(result)
    }
    result=sum(input[(row-radius):(row+radius),(col-radius):(col+radius)])-lq.calc(input[(row-radius):(row+radius),(col-radius):(col+radius)])*(2*radius+1)^2
    return(result)
  }
  cart2pol <- function(data,centroid){
    data.pol=data.frame('r'=sqrt((data$x-centroid[1])^2+(data$y-centroid[2])^2),'theta'=atan2((data$y-centroid[2]),(data$x-centroid[1]))*180/pi)
    return(data.pol)
  }
  pol2cart <- function(data,centroid){
    data.cart=data.frame('x'=data$r*cos(data$theta*pi/180)+centroid[1],'y'=data$r*sin(data$theta*pi/180)+centroid[2])
    return(data.cart)
  }
  cart2mat <- function(data,row.num){
    data.mat=data.frame('row'=round(row.num-data$y+1),'col'=round(data$x))
    return(data.mat)
  }
  mat2cart <- function(data,row.num){
    data.cart=data.frame('x'=data$col,'y'=row.num-data$row+1)
    return(data.cart)
  }
  pol.adjust <- function(input,r.factor,theta.factor){
    data.adj=data.frame('r'=data$r*r.factor,'theta'=data$theta+theta.factor)
    return(data.adj)
  }
  C5toC3 <- function(input,centroid,row.num,r.factor,theta.factor){
    cart2pol <- function(data,centroid){
      data.pol=data.frame('r'=sqrt((data$x-centroid[1])^2+(data$y-centroid[2])^2),'theta'=atan2((data$y-centroid[2]),(data$x-centroid[1]))*180/pi)
      return(data.pol)
    }
    pol2cart <- function(data,centroid){
      data.cart=data.frame('x'=data$r*cos(data$theta*pi/180)+centroid[1],'y'=data$r*sin(data$theta*pi/180)+centroid[2])
      return(data.cart)
    }
    cart2mat <- function(data,row.num){
      data.mat=data.frame('row'=round(row.num-data$y+1),'col'=round(data$x))
      return(data.mat)
    }
    mat2cart <- function(data,row.num){
      data.cart=data.frame('x'=data$col,'y'=row.num-data$row+1)
      return(data.cart)
    }
    pol.adjust <- function(data,r.factor,theta.factor){
      data.adj=data.frame('r'=data$r*r.factor,'theta'=data$theta+theta.factor)
      return(data.adj)
    }
    data=input
    data[,1:2]=pol2cart(pol.adjust(cart2pol(input[,1:2],centroid),r.factor,theta.factor),centroid)
    data[,3:4]=cart2mat(data[,1:2],row.num)
    return(data)
  }
  C3toC5 <- function(input,centroid,row.num,r.factor,theta.factor){
    cart2pol <- function(data,centroid){
      data.pol=data.frame('r'=sqrt((data$x-centroid[1])^2+(data$y-centroid[2])^2),'theta'=atan2((data$y-centroid[2]),(data$x-centroid[1]))*180/pi)
      return(data.pol)
    }
    pol2cart <- function(data,centroid){
      data.cart=data.frame('x'=data$r*cos(data$theta*pi/180)+centroid[1],'y'=data$r*sin(data$theta*pi/180)+centroid[2])
      return(data.cart)
    }
    cart2mat <- function(data,row.num){
      data.mat=data.frame('row'=round(row.num-data$y+1),'col'=round(data$x))
      return(data.mat)
    }
    mat2cart <- function(data,row.num){
      data.cart=data.frame('x'=data$col,'y'=row.num-data$row+1)
      return(data.cart)
    }
    pol.adjust <- function(data,r.factor,theta.factor){
      data.adj=data.frame('r'=data$r*r.factor,'theta'=data$theta+theta.factor)
      return(data.adj)
    }
    data=input
    data[,1:2]=pol2cart(pol.adjust(cart2pol(input[,1:2],centroid),1/r.factor,0-theta.factor),centroid)
    data[,3:4]=cart2mat(data[,1:2],row.num)
    return(data)
  }
  pair.finder <- function(ref,data,delta){
    x=data[,1]-ref[1]
    y=data[,2]-ref[2]
    r=sqrt(x^2+y^2)
    if (min(r)<=delta){
      id=which.min(r)
      return(id)
    }
    if (min(r)>delta){
      return(NULL)
    }
  }
  find.spots <- function(data,box.size,low.lim,high.lim,fill.radius){
    lq.calc <- function(input){
      input.2=na.omit(c(input))
      result=(input.2[order(input.2)])[round(0.25*length(input.2))]
      return(result)
    }
    vol.calc <- function(input,x,y,radius){
      lq.calc <- function(input){
        input.2=na.omit(c(input))
        result=(input.2[order(input.2)])[round(0.25*length(input.2))]
        return(result)
      }
      result=sum(input[(x-radius):(x+radius),(y-radius):(y+radius)])-lq.calc(input[(x-radius):(x+radius),(y-radius):(y+radius)])*(2*radius+1)^2
      return(result)
    }
    col.id=t(array(1:dim(data)[2],dim = dim(data)[c(2,1)]))
    row.id=array(1:dim(data)[1],dim = dim(data))
    peeks=array(F,dim = dim(col.id))
    volumes=array(0,dim = dim(col.id))
    for (x in (1+max(c(box.size,fill.radius))):(dim(peeks)[1]-max(c(box.size,fill.radius)))){
      for (y in (1+max(c(box.size,fill.radius))):(dim(peeks)[2]-max(c(box.size,fill.radius)))){
        peeks[x,y]=(data[x,y]==max(data[(x-box.size):(x+box.size),(y-box.size):(y+box.size)]) & sum(data[(x-box.size):(x+box.size),(y-box.size):(y+box.size)]==max(data[(x-box.size):(x+box.size),(y-box.size):(y+box.size)]))==1)
        if (peeks[x,y]==TRUE){
          volumes[x,y]=vol.calc(data,x,y,fill.radius)
        }
      }
    }
    rows=c(row.id[peeks==T])
    cols=c(col.id[peeks==T])
    vols=c(volumes[peeks==T])
    spots=data.frame('x'=cols[(vols>=low.lim & vols<=high.lim)],'y'=((dim(data)[1]+1)-rows)[(vols>=low.lim & vols<=high.lim)],'row'=rows[(vols>=low.lim & vols<=high.lim)],'col'=cols[(vols>=low.lim & vols<=high.lim)],'vol'=vols[(vols>=low.lim & vols<=high.lim)])
    final=spots[order(spots$vol),]
    return(final)
  }
  #
  if (is.null(Cy3.file.name)==TRUE){
    Cy3.file.name=list.files(path = path.to.file,pattern = '[.]tif')[1]
    show(paste0('Cy3 Filename = ',Cy3.file.name))
  }
  if (is.null(Cy5.file.name)==TRUE){
    Cy5.file.name=list.files(path = path.to.file,pattern = '[.]tif')[2]
    show(paste0('Cy5 Filename = ',Cy5.file.name))
  }
  Cy3.image.raw=suppressWarnings(tiff::readTIFF(paste0(path.to.file,Cy3.file.name),all = TRUE,as.is=TRUE))
  Cy5.image.raw=suppressWarnings(tiff::readTIFF(paste0(path.to.file,Cy5.file.name),all = TRUE,as.is=TRUE))
  pixel.size=dim(Cy3.image.raw[[1]])
  frame.number=length(Cy3.image.raw)
  Cy3.image.data=array(NA,dim = c(pixel.size[1],pixel.size[2],frame.number)); for (i in 1:frame.number){Cy3.image.data[,,i]=Cy3.image.raw[[i]]}
  Cy5.image.data=array(NA,dim = c(pixel.size[1],pixel.size[2],frame.number)); for (i in 1:frame.number){Cy5.image.data[,,i]=Cy5.image.raw[[i]]}
  rm(Cy3.image.raw,Cy5.image.raw)
  Cy3.image.avg=apply(Cy3.image.data,MARGIN = c(1,2),FUN = mean)
  Cy5.image.avg=apply(Cy5.image.data,MARGIN = c(1,2),FUN = mean)
  #
  load(file = paste0(path.to.file,align.file.name))
  r.factor=r.adj.value
  theta.factor=theta.adj.value
  row.num=dim(Cy3.image.avg)[1]
  centroid=rev(dim(Cy3.image.avg)/2+0.5)
  #
  par(fig = c(0,0.5,0,1))
  suppressWarnings(image(t(pracma::flipud(Cy3.image.avg)),col=gray.colors(cumprod(pixel.size)),x = 1:pixel.size[2],y=1:pixel.size[1],axes=FALSE,xlab='',ylab='',main = 'Cy3'))
  par(fig = c(0.5,1,0,1),new  = TRUE)
  suppressWarnings(image(t(pracma::flipud(Cy5.image.avg)),col=gray.colors(cumprod(pixel.size)),x = 1:pixel.size[2],y=1:pixel.size[1],axes=FALSE,xlab='',ylab='',main = 'Cy5'))
  check.1=readline(prompt = 'Image loading complete. Proceed with spot detection? (y/n):  '); if (check.1=='n'){stop('SCRIPT ABORTED BY USER')}; if (check.1!='n' & check.1!='y'){stop('INVALID RESPONSE')}
  #
  Cy3.spots=find.spots(data = Cy3.image.avg,box.size = Cy3.search.box,fill.radius = Cy3.integration.radius,low.lim = Cy3.spot.min,high.lim = Cy3.spot.max)
  Cy5.spots=find.spots(data = Cy5.image.avg,box.size = Cy5.search.box,fill.radius = Cy5.integration.radius,low.lim = Cy5.spot.min,high.lim = Cy5.spot.max)
  Cy5spots.Cy5axes=Cy5.spots
  Cy5spots.Cy3axes=C5toC3(input = Cy5.spots,centroid = centroid,row.num = row.num,r.factor = r.factor,theta.factor = theta.factor)
  Cy3spots.Cy3axes=Cy3.spots
  Cy3spots.Cy5axes=C3toC5(input = Cy3.spots,centroid = centroid,row.num = row.num,r.factor = r.factor,theta.factor = theta.factor)
  par(fig = c(0,0.5,0,1))
  suppressWarnings(image(t(pracma::flipud(Cy3.image.avg)),col=gray.colors(cumprod(pixel.size)),x = 1:pixel.size[2],y=1:pixel.size[1],axes=FALSE,xlab='',ylab='',main = 'Cy3'))
  points(Cy5spots.Cy3axes$x,Cy5spots.Cy3axes$y,col='red',cex=0.5,pch='x')
  points(Cy3.spots$x,Cy3.spots$y,col='green',cex=1.2)
  show(paste0('Median Cy3 spot intensity = ',median(Cy3.spots$vol)))
  show(paste0('Number of Cy3 spots identified = ',length(Cy3.spots$x)))
  par(fig = c(0.5,1,0,1),new  = TRUE)
  suppressWarnings(image(t(pracma::flipud(Cy5.image.avg)),col=gray.colors(cumprod(pixel.size)),x = 1:pixel.size[2],y=1:pixel.size[1],axes=FALSE,xlab='',ylab='',main = 'Cy5'))
  points(Cy3spots.Cy5axes$x,Cy3spots.Cy5axes$y,col='green',cex=0.5,pch='x')
  points(Cy5.spots$x,Cy5.spots$y,col='red',cex=1.2)
  show(paste0('Median Cy5 spot intensity = ',median(Cy5.spots$vol)))
  show(paste0('Number of Cy5 spots identified = ',length(Cy5.spots$x)))
  check.2=readline(prompt = 'Is spot selection acceptable -- proceed to spot alignment and pairing? (y/n/quit):  '); if (check.2=='quit'){stop('SCRIPT ABORTED BY USER')}; if (check.2!='n' & check.2!='y' & check.2!='quit'){stop('INVALID RESPONSE')}; while (check.2=='n'){
    Cy3.search.box=as.numeric(readline(prompt = 'New Cy3 spot search box radius (ENTER for default):  '))
    if (is.na(Cy3.search.box)==T){
      Cy3.search.box=6
    }
    Cy5.search.box=as.numeric(readline(prompt = 'New Cy5 spot search box radius (ENTER for default):  '))
    if (is.na(Cy5.search.box)==T){
      Cy5.search.box=6
    }
    Cy3.integration.radius=as.numeric(readline(prompt = 'New Cy3 spot integration box radius (ENTER for default):  '))
    if (is.na(Cy3.integration.radius)==T){
      Cy3.integration.radius=6
    }
    Cy5.integration.radius=as.numeric(readline(prompt = 'New Cy5 spot integration box radius (ENTER for default):  '))
    if (is.na(Cy5.integration.radius)==T){
      Cy5.integration.radius=6
    }
    Cy3.spot.min=as.numeric(readline(prompt = 'New Cy3 spot volume minimum (ENTER for default):  '))
    if (is.na(Cy3.spot.min)==T){
      Cy3.spot.min=3500
    }
    Cy5.spot.min=as.numeric(readline(prompt = 'New Cy5 spot volume minimum (ENTER for default):  '))
    if (is.na(Cy5.spot.min)==T){
      Cy5.spot.min=1e4
    }
    Cy3.spot.max=as.numeric(readline(prompt = 'New Cy3 spot volume maximum (ENTER for default):  '))
    if (is.na(Cy3.spot.max)==T){
      Cy3.spot.max=Inf
    }
    Cy5.spot.max=as.numeric(readline(prompt = 'New Cy5 spot volume maximum (ENTER for default):  '))
    if (is.na(Cy5.spot.max)==T){
      Cy5.spot.max=Inf
    }
    Cy3.spots=find.spots(data = Cy3.image.avg,box.size = Cy3.search.box,fill.radius = Cy3.integration.radius,low.lim = Cy3.spot.min,high.lim = Cy3.spot.max)
    Cy5.spots=find.spots(data = Cy5.image.avg,box.size = Cy5.search.box,fill.radius = Cy5.integration.radius,low.lim = Cy5.spot.min,high.lim = Cy5.spot.max)
    Cy5spots.Cy5axes=Cy5.spots
    Cy5spots.Cy3axes=C5toC3(input = Cy5.spots,centroid = centroid,row.num = row.num,r.factor = r.factor,theta.factor = theta.factor)
    Cy3spots.Cy3axes=Cy3.spots
    Cy3spots.Cy5axes=C3toC5(input = Cy3.spots,centroid = centroid,row.num = row.num,r.factor = r.factor,theta.factor = theta.factor)
    par(fig = c(0,0.5,0,1))
    suppressWarnings(image(t(pracma::flipud(Cy3.image.avg)),col=gray.colors(cumprod(pixel.size)),x = 1:pixel.size[2],y=1:pixel.size[1],axes=FALSE,xlab='',ylab='',main = 'Cy3'))
    points(Cy5spots.Cy3axes$x,Cy5spots.Cy3axes$y,col='red',cex=0.5,pch='x')
    points(Cy3.spots$x,Cy3.spots$y,col='green',cex=1.2)
    show(paste0('Median Cy3 spot intensity = ',median(Cy3.spots$vol)))
    show(paste0('Number of Cy3 spots identified = ',length(Cy3.spots$x)))
    par(fig = c(0.5,1,0,1),new  = TRUE)
    suppressWarnings(image(t(pracma::flipud(Cy5.image.avg)),col=gray.colors(cumprod(pixel.size)),x = 1:pixel.size[2],y=1:pixel.size[1],axes=FALSE,xlab='',ylab='',main = 'Cy5'))
    points(Cy3spots.Cy5axes$x,Cy3spots.Cy5axes$y,col='green',cex=0.5,pch='x')
    points(Cy5.spots$x,Cy5.spots$y,col='red',cex=1.2)
    show(paste0('Median Cy5 spot intensity = ',median(Cy5.spots$vol)))
    show(paste0('Number of Cy5 spots identified = ',length(Cy5.spots$x)))
    check.2=readline(prompt = 'Is spot selection acceptable -- proceed to spot alignment and pairing? (y/n/quit):  ')
    if (check.2=='quit'){stop('SCRIPT ABORTED BY USER')}
    if (check.2!='n' & check.2!='y' & check.2!='quit'){stop('INVALID RESPONSE')}
  }
  #
  delta=(0:20)
  Cy3.pairs=matrix(0,nrow=nrow(Cy3spots.Cy3axes),ncol = length(delta))
  Cy5.pairs=matrix(0,nrow=nrow(Cy5spots.Cy3axes),ncol = length(delta))
  for (i in 1:length(delta)){
    COUNTER=0
    for (j in 1:nrow(Cy5spots.Cy3axes)){
      pair.id=pair.finder(ref = Cy5spots.Cy3axes[j,1:2],data = Cy3spots.Cy3axes[,1:2],delta = delta[i])
      if (is.null(pair.id)==FALSE){
        COUNTER=COUNTER+1
        Cy5.pairs[j,i]=COUNTER
        Cy3.pairs[pair.id,i]=COUNTER
      }
    }
  }
  rm(COUNTER)
  #
  delta.info=data.frame('delta'=delta,'Cy3.Spots'=colSums(Cy3.pairs==0),'Cy5.Spots'=colSums(Cy5.pairs==0),'DUAL.Spots'=colSums(Cy5.pairs!=0))
  show(delta.info)
  check.3='n'; while (check.3=='n'){
    input.1=as.numeric(readline(prompt = 'Which delta value would you like to proceed with? (0-10,quit):  '))
    temp.1.Cy3axes=rbind(Cy3spots.Cy3axes[Cy3.pairs[,delta==input.1]==0,1:4],Cy5spots.Cy3axes[Cy5.pairs[,delta==input.1]==0,1:4],(Cy3spots.Cy3axes[((Cy3.pairs[,delta==input.1])!=0),1:4])[order(Cy3.pairs[((Cy3.pairs[,delta==input.1])!=0),delta==input.1]),])
    temp.1.Cy5axes=rbind(Cy3spots.Cy5axes[Cy3.pairs[,delta==input.1]==0,1:4],Cy5spots.Cy5axes[Cy5.pairs[,delta==input.1]==0,1:4],Cy5spots.Cy5axes[Cy5.pairs[,delta==input.1]!=0,1:4])
    ALLspots.Cy3axes=data.frame('Particle.ID'=1:nrow(temp.1.Cy3axes),'x'=temp.1.Cy3axes[,1],'y'=temp.1.Cy3axes[,2],'row'=temp.1.Cy3axes[,3],'col'=temp.1.Cy3axes[,4],'Type.ID'=rep(c('Cy3','Cy5','DUAL'),times = delta.info[delta.info$delta==input.1,2:4]))
    ALLspots.Cy5axes=data.frame('Particle.ID'=1:nrow(temp.1.Cy5axes),'x'=temp.1.Cy5axes[,1],'y'=temp.1.Cy5axes[,2],'row'=temp.1.Cy5axes[,3],'col'=temp.1.Cy5axes[,4],'Type.ID'=rep(c('Cy3','Cy5','DUAL'),times = delta.info[delta.info$delta==input.1,2:4]))
    par(fig = c(0,0.5,0,1))
    suppressWarnings(image(t(pracma::flipud(Cy3.image.avg)),col=gray.colors(cumprod(pixel.size)),x = 1:pixel.size[2],y=1:pixel.size[1],axes=FALSE,xlab='',ylab='',main = 'Cy3'))
    points(ALLspots.Cy3axes$x,ALLspots.Cy3axes$y,col=rep(c('green','red','blue'),times = delta.info[delta.info$delta==input.1,2:4]),cex=1.2)
    par(fig = c(0.5,1,0,1),new  = TRUE)
    suppressWarnings(image(t(pracma::flipud(Cy5.image.avg)),col=gray.colors(cumprod(pixel.size)),x = 1:pixel.size[2],y=1:pixel.size[1],axes=FALSE,xlab='',ylab='',main = 'Cy5'))
    points(ALLspots.Cy5axes$x,ALLspots.Cy5axes$y,col=rep(c('green','red','blue'),times = delta.info[delta.info$delta==input.1,2:4]),cex=1.2)
    check.3=readline(prompt = 'Is delta selection acceptable -- proceed to calculations and file export? (y/n/quit):  ')
    if (check.3=='quit'){stop('SCRIPT ABORTED BY USER')}
    if (check.3!='n' & check.3!='y' & check.3!='quit'){stop('INVALID RESPONSE')}
  }
  #
  Cy3.traces=matrix(0,nrow = frame.number,ncol = nrow(ALLspots.Cy3axes))
  Cy5.traces=matrix(0,nrow = frame.number,ncol = nrow(ALLspots.Cy5axes))
  Cy3.snaps=array(0,dim = c(2*Cy3.integration.radius+1,2*Cy3.integration.radius+1,10,nrow(ALLspots.Cy3axes)))
  Cy5.snaps=array(0,dim = c(2*Cy5.integration.radius+1,2*Cy5.integration.radius+1,10,nrow(ALLspots.Cy5axes)))
  snap.length=frame.number/10
  for (i in 1:nrow(ALLspots.Cy3axes)){
    Cy3.traces[,i]=apply(Cy3.image.data,MARGIN = c(3),FUN = vol.calc,row=ALLspots.Cy3axes$row[i],col=ALLspots.Cy3axes$col[i],radius=Cy3.integration.radius)
    Cy5.traces[,i]=apply(Cy5.image.data,MARGIN = c(3),FUN = vol.calc,row=ALLspots.Cy5axes$row[i],col=ALLspots.Cy5axes$col[i],radius=Cy5.integration.radius)
    for (j in 1:10){
      Cy3.snaps[,,j,i]=apply(Cy3.image.data[(ALLspots.Cy3axes$row[i]-Cy3.integration.radius):(ALLspots.Cy3axes$row[i]+Cy3.integration.radius),(ALLspots.Cy3axes$col[i]-Cy3.integration.radius):(ALLspots.Cy3axes$col[i]+Cy3.integration.radius),((j-1)*snap.length+1):(j*snap.length)],MARGIN = c(1,2),FUN = mean)
      Cy5.snaps[,,j,i]=apply(Cy5.image.data[(ALLspots.Cy5axes$row[i]-Cy5.integration.radius):(ALLspots.Cy5axes$row[i]+Cy5.integration.radius),(ALLspots.Cy5axes$col[i]-Cy5.integration.radius):(ALLspots.Cy5axes$col[i]+Cy5.integration.radius),((j-1)*snap.length+1):(j*snap.length)],MARGIN = c(1,2),FUN = mean)
    }
  }
  row.names(Cy3.traces)=time.step*(1:frame.number)
  row.names(Cy5.traces)=time.step*(1:frame.number)
  id.spots.output.files=list('Cy3'=list('Cy3 composite image data'=Cy3.image.avg,'initial particle Cy3-summary'=ALLspots.Cy3axes,'Cy3 traces'=Cy3.traces,'Cy3 particle images data'=Cy3.snaps),'Cy5'=list('Cy5 composite image data'=Cy5.image.avg,'initial particle Cy5-summary'=ALLspots.Cy5axes,'Cy5 traces'=Cy5.traces,'Cy5 particle images data'=Cy5.snaps))
  #
  dev.print(pdf,paste0(path.to.file,'Identified_Spots.pdf'))
  utils::write.table(ALLspots.Cy3axes,file = paste0(path.to.file,'initial_particle_Cy3-summary.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  utils::write.table(ALLspots.Cy5axes,file = paste0(path.to.file,'initial_particle_Cy5-summary.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  utils::write.table(Cy3.traces,file = paste0(path.to.file,'initial_particle_Cy3-traces.txt'),quote = FALSE,sep = '\t',row.names = TRUE,col.names = FALSE)
  utils::write.table(Cy5.traces,file = paste0(path.to.file,'initial_particle_Cy5-traces.txt'),quote = FALSE,sep = '\t',row.names = TRUE,col.names = FALSE)
  save(list = c('Cy3.image.avg','ALLspots.Cy3axes','Cy3.traces','Cy3.snaps','Cy3.integration.radius','Cy5.image.avg','ALLspots.Cy5axes','Cy5.traces','Cy5.snaps','Cy5.integration.radius','time.step','frame.number','pixel.size'),file = paste0(path.to.file,'Initial-Particle_Data.RData'))
  return(id.spots.output.files)
}

FRET.refine <- function(path.to.file='./',file.name='Initial-Particle_Data.RData',auto.filter='none',spot.types='all',classification.strategy='classic',background.subtraction='lower.quartile'){
  count.condense <- function(input){
    new.count=input
    for (i in 1:length(input)){
      new.count[i]=length(unique(input[1:i]))
    }
    return(new.count)
  }
  woh.scale <- function(data){
    a=c(data)
    b=(data-min(a))/diff(range(a))
    return(b)
  }
  lq.calc <- function(input){
    input.2=na.omit(c(input))
    result=(input.2[order(input.2)])[round(0.25*length(input.2))]
    return(result)
  }
  k.mode <- function(data){
    k=density(data)
    mode=k[['x']][which.max(k[['y']])]
    return(mode)
  }
  classify.states <- function(data,signal.step,fun='classic',background=0){
    if (fun=='classic'){
      fit=smooth.spline((1:length(data))[is.na(data)==F],na.omit(data),spar = 0.5)
      d1=predict(object = fit,x = 1:length(data),deriv = 1)
      swaps=c(1,d1[['x']][(abs(d1[['y']])>=1.0) & (splus2R::peaks(x = abs(d1[['y']]),span = 5))],max(d1[['x']]))
      pos.averages=rep(NA,times=length(data))
      states=rep(0,times=length(data))
      state.averages=rep(NA,times=length(data))
      for (i in 1:(length(swaps)-1)){
        pos.averages[swaps[i]:swaps[i+1]]=mean(na.omit(data[swaps[i]:swaps[i+1]]))
      }
      for (j in 0:50){
        states[(pos.averages>=(signal.step*j-signal.step/2+background)) & (pos.averages<=(signal.step*j+signal.step/2+background))]=j
      }
      for (k in 0:50){
        state.averages[states==k]=mean(na.omit(data[states==k]))
      }
      final=as.data.frame(cbind(states,state.averages))
      colnames(final)=c('id','avg')
      return(final)
    }
    if (fun=='basic'){
      pos.averages=matrix(NA,ncol = ncol(data),nrow = nrow(data))
      for (n in 1:ncol(data)){
        fit=smooth.spline((1:length(data[,n]))[is.na(data[,n])==F],na.omit(data[,n]),spar = 0.5)
        d1=predict(object = fit,x = 1:length(data[,n]),deriv = 1)
        swaps=c(1,d1[['x']][(abs(d1[['y']])>=1.0) & (splus2R::peaks(x = abs(d1[['y']]),span = 5))],max(d1[['x']]))
        for (i in 1:(length(swaps)-1)){
          pos.averages[swaps[i]:swaps[i+1],n]=mean(na.omit(data[swaps[i]:swaps[i+1],n]))
        }
      }
      final=pos.averages
      return(final)
    }
    if (fun=='uni.dbscan'){
      pos.averages=matrix(NA,ncol = ncol(data),nrow = nrow(data))
      for (n in 1:ncol(data)){
        fit=smooth.spline((1:length(data[,n]))[is.na(data[,n])==F],na.omit(data[,n]),spar = 0.5)
        d1=predict(object = fit,x = 1:length(data[,n]),deriv = 1)
        swaps=c(1,d1[['x']][(abs(d1[['y']])>=1.0) & (splus2R::peaks(x = abs(d1[['y']]),span = 5))],max(d1[['x']]))
        for (i in 1:(length(swaps)-1)){
          pos.averages[swaps[i]:swaps[i+1],n]=mean(na.omit(data[swaps[i]:swaps[i+1],n]))
        }
      }
      uni.dbscan <- function(datter,eps){
        temp.0=order(c(datter))
        temp.1=datter[temp.0]
        temp.2=c(abs(temp.1[1:(length(temp.1)-1)]-temp.1[2:length(temp.1)]),0)>eps
        temp.3=c(1,(1:length(temp.2))[temp.2],length(temp.2))
        temp.4=rep(NA,times=length(temp.2))
        for (i in 0:(length(temp.3)-2)){
          temp.4[(temp.3[i+1]):(temp.3[i+2])]=i
        }
        temp.5=matrix(temp.4[order(temp.0)],nrow = nrow(datter),ncol = ncol(datter))
        return(temp.5)
      }
      classes=uni.dbscan(pos.averages,eps = 1*sd(na.omit(data-pos.averages)))
      values=matrix(NA,nrow = nrow(classes),ncol = ncol(classes))
      for (k in 0:max(classes)){
        values[classes==k]=mean(na.omit(data[classes==k]))
      }
      final=list('id'=classes,'avg'=values)
      return(final)
    }
    if (fun=='var.shift'){
      pos.averages=rep(NA,times=length(data))
      fit=smooth.spline((1:length(data))[is.na(data)==F],na.omit(data),spar = 0.5)
      d1=predict(object = fit,x = 1:length(data),deriv = 1)
      swaps=c(1,d1[['x']][(abs(d1[['y']])>=1.0) & (splus2R::peaks(x = abs(d1[['y']]),span = 5))],max(d1[['x']]))
      for (i in 1:(length(swaps)-1)){
        pos.averages[swaps[i]:swaps[i+1]]=mean(na.omit(data[swaps[i]:swaps[i+1]]))
      }
      var.shift <- function(datter,eps){
        temp.0=order(datter)
        temp.1=datter[temp.0]
        temp.2=c(abs(temp.1[1:(length(temp.1)-1)]-temp.1[2:length(temp.1)]),0)>eps
        temp.3=c(1,(1:length(temp.2))[temp.2],length(temp.2))
        temp.4=rep(NA,times=length(temp.2))
        for (i in 1:(length(temp.3)-1)){
          temp.4[(temp.3[i]):(temp.3[i+1])]=i
        }
        temp.5=temp.4[order(temp.0)]
        return(temp.5)
      }
      classes=var.shift(pos.averages,eps = 1.7*sd(na.omit(data-pos.averages)))
      values=rep(NA,times=length(classes))
      for (k in 0:max(classes)){
        values[classes==k]=mean(na.omit(data[classes==k]))
      }
      classes[abs(values)<=2*sd(na.omit(data-pos.averages))]=0
      final=data.frame('id'=classes,'avg'=values)
      return(final)
    }
  }
  #
  load(file = paste0(path.to.file,file.name))
  Cy3.trace.rolls=apply(Cy3.traces,MARGIN = c(2),FUN = data.table::frollmean,n=5,align='center')
  Cy5.trace.rolls=apply(Cy5.traces,MARGIN = c(2),FUN = data.table::frollmean,n=5,align='center')
  if (background.subtraction=='k.mode'){
    Cy3.background=k.mode(na.omit(Cy3.trace.rolls))
    Cy3.trace.rolls=Cy3.trace.rolls-Cy3.background
    Cy5.background=k.mode(na.omit(Cy5.trace.rolls))
    Cy5.trace.rolls=Cy5.trace.rolls-Cy5.background
  }
  if (background.subtraction=='basal.states'){
    Cy3.background=median(na.omit(apply(classify.states(Cy3.trace.rolls,fun = 'basic'),MARGIN = c(2),FUN = min)))
    Cy3.trace.rolls=Cy3.trace.rolls-Cy3.background
    Cy5.background=median(na.omit(apply(classify.states(Cy5.trace.rolls,fun = 'basic'),MARGIN = c(2),FUN = min)))
    Cy5.trace.rolls=Cy5.trace.rolls-Cy5.background
  }
  if (background.subtraction=='lower.quartile'){
    Cy3.background=lq.calc(classify.states(Cy3.trace.rolls,fun = 'basic'))
    Cy3.trace.rolls=Cy3.trace.rolls-Cy3.background
    Cy5.background=lq.calc(classify.states(Cy5.trace.rolls,fun = 'basic'))
    Cy5.trace.rolls=Cy5.trace.rolls-Cy5.background
  }
  #
  particles.to.keep=rep(FALSE,times=nrow(ALLspots.Cy3axes))
  residence.times=c('Particle','Start','Stop','Residence')
  Cy3.dwell.calls=c('Particle','State','Dwell')
  Cy3.residence.calls=c('Particle','Bound','Residence')
  Cy5.dwell.calls=c('Particle','State','Dwell')
  Cy5.residence.calls=c('Particle','Bound','Residence')
  if (classification.strategy=='classic' | classification.strategy=='var.shift'){
    Cy3.state.calls=matrix(NA,nrow = frame.number,ncol = nrow(ALLspots.Cy3axes))
    Cy5.state.calls=matrix(NA,nrow = frame.number,ncol = nrow(ALLspots.Cy5axes))
    Cy3.signal.step=3*sd(na.omit(Cy3.trace.rolls-classify.states(Cy3.trace.rolls,fun = 'basic')))
    Cy5.signal.step=3*sd(na.omit(Cy5.trace.rolls-classify.states(Cy5.trace.rolls,fun = 'basic')))
  }
  if (classification.strategy=='uni.dbscan'){
    Cy3.classifications=classify.states(Cy3.trace.rolls,fun = 'uni.dbscan')
    Cy3.state.calls=Cy3.classifications[['id']]
    Cy5.classifications=classify.states(Cy5.trace.rolls,fun = 'uni.dbscan')
    Cy5.state.calls=Cy5.classifications[['id']]
  }
  Cy3.filtering=rep(F,times=nrow(ALLspots.Cy3axes))
  Cy5.filtering=rep(F,times=nrow(ALLspots.Cy5axes))
  for (i in 1:nrow(ALLspots.Cy3axes)){
    if (classification.strategy=='classic' | classification.strategy=='var.shift'){
      Cy3.states=classify.states(Cy3.trace.rolls[,i],Cy3.signal.step,fun = classification.strategy)
      Cy3.state.calls[,i]=Cy3.states$id
      Cy5.states=classify.states(Cy5.trace.rolls[,i],Cy5.signal.step,fun = classification.strategy)
      Cy5.state.calls[,i]=Cy5.states$id
    }
    Cy3.dwell.calls=rbind(Cy3.dwell.calls,cbind(rep(i,times=length(rle(Cy3.state.calls[,i])[['values']])),rle(Cy3.state.calls[,i])[['values']],rle(Cy3.state.calls[,i])[['lengths']]*time.step))
    Cy3.residence.calls=rbind(Cy3.residence.calls,cbind(rep(i,times=length(rle(Cy3.state.calls[,i]>=1)[['values']])),rle(Cy3.state.calls[,i]>=1)[['values']],rle(Cy3.state.calls[,i]>=1)[['lengths']]*time.step))
    Cy5.dwell.calls=rbind(Cy5.dwell.calls,cbind(rep(i,times=length(rle(Cy5.state.calls[,i])[['values']])),rle(Cy5.state.calls[,i])[['values']],rle(Cy5.state.calls[,i])[['lengths']]*time.step))
    Cy5.residence.calls=rbind(Cy5.residence.calls,cbind(rep(i,times=length(rle(Cy5.state.calls[,i]>=1)[['values']])),rle(Cy5.state.calls[,i]>=1)[['values']],rle(Cy5.state.calls[,i]>=1)[['lengths']]*time.step))
    if  (auto.filter=='unbound') {
      Cy3.filtering[i]=(sum(Cy3.state.calls[,i])==0)
      Cy5.filtering[i]=(sum(Cy5.state.calls[,i])==0)
    }
    if  (auto.filter=='all.stable') {
      Cy3.filtering[i]=(length(table(Cy3.state.calls[,i]))==1)
      Cy5.filtering[i]=(length(table(Cy5.state.calls[,i]))==1)
    }
  }
  #
  par(fig=c(0,1,0,1))
  plot((1:frame.number)*time.step,apply(Cy3.state.calls>0,MARGIN = c(1),FUN = sum)/nrow(ALLspots.Cy3axes),type='l',col='green',main='Photobleaching Check',xlab='Time (s)',ylab='Proportion of Bound Particles',ylim = range(c(apply(Cy3.state.calls>0,MARGIN = c(1),FUN = sum)/nrow(ALLspots.Cy3axes),apply(Cy5.state.calls>0,MARGIN = c(1),FUN = sum)/nrow(ALLspots.Cy5axes))))
  lines((1:frame.number)*time.step,apply(Cy5.state.calls>0,MARGIN = c(1),FUN = sum)/nrow(ALLspots.Cy5axes),col='red')
  check.1=readline('Proceed with particle refinement? (y/n):  ')
  if (check.1=='n'){stop('SCRIPT ABORTED BY USER')}
  if (check.1!='n' & check.1!='y'){stop('INVALID RESPONSE')}
  #
  COUNTER=0
  if (spot.types=='all'){
    desired.type=rep(TRUE,times=nrow(ALLspots.Cy3axes))
  } else {
    desired.type=(ALLspots.Cy3axes$Type.ID==spot.types)
  }
  for (i in 1:nrow(ALLspots.Cy3axes)){
    if ((Cy3.filtering[i]*Cy3.filtering[i])==FALSE & desired.type[i]==TRUE){
      if (classification.strategy=='classic' | classification.strategy=='var.shift'){
        Cy3.states=classify.states(Cy3.trace.rolls[,i],Cy3.signal.step,fun = classification.strategy)
        Cy5.states=classify.states(Cy5.trace.rolls[,i],Cy5.signal.step,fun = classification.strategy)
      }
      if (classification.strategy=='uni.dbscan'){
        Cy3.states=data.frame('id'=Cy3.classifications[['id']][,i],'avg'=Cy3.classifications[['avg']][,i])
        Cy5.states=data.frame('id'=Cy5.classifications[['id']][,i],'avg'=Cy5.classifications[['avg']][,i])
      }
      par(fig=c(0,1,0.7,1))
      temp.1=matrix(c(woh.scale(Cy3.snaps[,,,i]),woh.scale(Cy5.snaps[,,,i])),nrow = (2*Cy3.integration.radius+1),ncol = (2*Cy3.integration.radius+1)*20)
      image(t(pracma::flipud(rbind(temp.1[,1:(ncol(temp.1)/2)],temp.1[,(ncol(temp.1)/2+1):ncol(temp.1)]))),col=gray.colors(length(temp.1)),axes=FALSE)
      par(fig=c(0,1,0.5,0.75),new = TRUE)
      plot((1:frame.number)*time.step,Cy3.trace.rolls[,i],type='l',main = 'Cy3 Trace',xlab='Time (s)',ylab = 'Signal')
      lines(((1:frame.number)*time.step)[is.na(Cy3.states$avg)==F],Cy3.states$avg[is.na(Cy3.states$avg)==F],col='green',lwd=3)
      par(fig=c(0,1,0.25,0.5),new = TRUE)
      plot((1:frame.number)*time.step,Cy5.trace.rolls[,i],type='l',main = 'Cy5 Trace',xlab='Time (s)',ylab = 'Signal')
      lines(((1:frame.number)*time.step)[is.na(Cy5.states$avg)==F],Cy5.states$avg[is.na(Cy5.states$avg)==F],col='red',lwd=3)
      par(fig=c(0,0.5,0,0.25),new = TRUE)
      plot(((1:frame.number)*time.step)[is.na(Cy3.states$avg)==F],Cy3.states$avg[is.na(Cy3.states$avg)==F],type='l',main = 'State Overlay',xlab='Time (s)',ylab = 'Signal',col='green',ylim = range(na.omit(c(Cy3.states$avg[is.na(Cy3.states$avg)==F],Cy5.states$avg[is.na(Cy5.states$avg)==F]))))
      lines(((1:frame.number)*time.step)[is.na(Cy5.states$avg)==F],Cy5.states$avg[is.na(Cy5.states$avg)==F],col='red')
      par(fig=c(0.5,1,0,0.25),new = TRUE)
      plot((1:frame.number)*time.step,Cy3.trace.rolls[,i],type='l',main = 'Trace Overlay',xlab='Time (s)',ylab = 'Signal',col='green',ylim = range(na.omit(c(Cy3.trace.rolls[,i],Cy5.trace.rolls[,i]))))
      lines((1:frame.number)*time.step,Cy5.trace.rolls[,i],col='red')
      check.2=readline('Should this particle be used for analysis? (y/n/quit/finish):  ')
      if (check.2=='y'){
        particles.to.keep[i]=TRUE
        COUNTER=COUNTER+1
        event.number=as.numeric(readline('How many events will you record for this particle?:  '))
        if (length(event.number)==1 & is.na(event.number)==FALSE){
          event.times=matrix(0,nrow=event.number,ncol = 4); event.times[,1]=rep(COUNTER,times=event.number)
          for (j in 1:event.number){
            show('Please click on the starting point then stopping point of a binding event (Trace Overlay -- x-axis).')
            event.times[j,2:3]=as.numeric(identify((1:frame.number)*time.step,rep(0,times=frame.number),n=2,plot = FALSE)*time.step)
            event.times[j,4]=diff(event.times[j,2:3])
            arrows(x0 = event.times[j,2],x1 = event.times[j,3],y0=max(na.omit(c(Cy3.trace.rolls[,i],Cy5.trace.rolls[,i]))),y1=max(na.omit(c(Cy3.trace.rolls[,i],Cy5.trace.rolls[,i]))),angle = 90,code = 3,col = 'black')
          }
          residence.times=rbind(residence.times,event.times)
        }
        dev.print(pdf,paste0(path.to.file,'Particle_Trace_',COUNTER,'.pdf'))
      }
      if (check.2=='quit'){
        stop('SCRIPT ABORTED BY USER')
      }
    }
    if (check.2=='finish'){
      break
    }
  }
  show(paste0('Particle refinement completed -- beginning data exports!'))
  if (length(dim(residence.times))==2){
    residence.data=as.data.frame(as.matrix(residence.times[2:nrow(residence.times),])); for (k in 1:4){residence.data[,k]=as.numeric(residence.data[,k])}; colnames(residence.data)=c('Particle','Start','Stop','Residence')
  } else {
    residence.data=c('Particle','Start','Stop','Residence')
  }
  Cy3.residence.calls=as.data.frame(as.matrix(Cy3.residence.calls[2:nrow(Cy3.residence.calls),])); for (k in 1:3){Cy3.residence.calls[,k]=as.numeric(Cy3.residence.calls[,k])}; colnames(Cy3.residence.calls)=c('Particle','Bound','Residence')
  Cy3.dwell.calls=as.data.frame(as.matrix(Cy3.dwell.calls[2:nrow(Cy3.dwell.calls),])); for (k in 1:3){Cy3.dwell.calls[,k]=as.numeric(Cy3.dwell.calls[,k])}; colnames(Cy3.dwell.calls)=c('Particle','State','Dwell')
  refined.Cy3.traces=Cy3.traces[,particles.to.keep]
  refined.Cy3.state.calls=Cy3.state.calls[,particles.to.keep]
  refined.Cy3.residence.calls=Cy3.residence.calls[which(Cy3.residence.calls[,1]%in%which(particles.to.keep)),]; refined.Cy3.residence.calls[,1]=count.condense(refined.Cy3.residence.calls[,1])
  refined.Cy3.dwell.calls=Cy3.dwell.calls[which(Cy3.dwell.calls[,1]%in%which(particles.to.keep)),]; refined.Cy3.dwell.calls[,1]=count.condense(refined.Cy3.dwell.calls[,1])
  refined.Cy3.trace.rolls=Cy3.trace.rolls[,particles.to.keep]
  Cy5.residence.calls=as.data.frame(as.matrix(Cy5.residence.calls[2:nrow(Cy5.residence.calls),])); for (k in 1:3){Cy5.residence.calls[,k]=as.numeric(Cy5.residence.calls[,k])}; colnames(Cy5.residence.calls)=c('Particle','Bound','Residence')
  Cy5.dwell.calls=as.data.frame(as.matrix(Cy5.dwell.calls[2:nrow(Cy5.dwell.calls),])); for (k in 1:3){Cy5.dwell.calls[,k]=as.numeric(Cy5.dwell.calls[,k])}; colnames(Cy5.dwell.calls)=c('Particle','State','Dwell')
  refined.Cy5.traces=Cy5.traces[,particles.to.keep]
  refined.Cy5.state.calls=Cy5.state.calls[,particles.to.keep]
  refined.Cy5.residence.calls=Cy5.residence.calls[which(Cy5.residence.calls[,1]%in%which(particles.to.keep)),]; refined.Cy5.residence.calls[,1]=count.condense(refined.Cy5.residence.calls[,1])
  refined.Cy5.dwell.calls=Cy5.dwell.calls[which(Cy5.dwell.calls[,1]%in%which(particles.to.keep)),]; refined.Cy5.dwell.calls[,1]=count.condense(refined.Cy5.dwell.calls[,1])
  refined.Cy5.trace.rolls=Cy5.trace.rolls[,particles.to.keep]
  REFspots.Cy3axes=ALLspots.Cy3axes[particles.to.keep,]
  refined.Cy3.snaps=Cy3.snaps[,,,particles.to.keep]
  REFspots.Cy5axes=ALLspots.Cy5axes[particles.to.keep,]
  refined.Cy5.snaps=Cy5.snaps[,,,particles.to.keep]
  #
  save(list = c('residence.data','Cy3.image.avg','refined.Cy3.traces','refined.Cy3.trace.rolls','REFspots.Cy3axes','refined.Cy3.snaps','Cy3.state.calls','Cy3.residence.calls','Cy3.dwell.calls','Cy5.image.avg','refined.Cy5.traces','refined.Cy5.trace.rolls','REFspots.Cy5axes','refined.Cy5.snaps','Cy5.state.calls','Cy5.residence.calls','Cy5.dwell.calls'),file = paste0(path.to.file,'Refined-Particle_Data.RData'))
  utils::write.table(residence.data,file = paste0(path.to.file,'residence_data.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  utils::write.table(refined.Cy3.traces,file = paste0(path.to.file,'selected_Cy3_traces.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  utils::write.table(refined.Cy3.trace.rolls,file = paste0(path.to.file,'selected_Cy3_smoothed-traces.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  utils::write.table(REFspots.Cy3axes,file = paste0(path.to.file,'selected_Cy3-particle_summary.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  utils::write.table(refined.Cy5.traces,file = paste0(path.to.file,'selected_Cy5_traces.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  utils::write.table(refined.Cy5.trace.rolls,file = paste0(path.to.file,'selected_Cy5_smoothed-traces.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  utils::write.table(REFspots.Cy5axes,file = paste0(path.to.file,'selected_Cy5-particle_summary.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  utils::write.table(Cy3.state.calls,file = paste0(path.to.file,'all-particle_Cy3-state-calls.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  utils::write.table(Cy3.dwell.calls,file = paste0(path.to.file,'all-particle_Cy3-dwell-calls.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  utils::write.table(Cy3.residence.calls,file = paste0(path.to.file,'all-particle_Cy3-residence-calls.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  utils::write.table(Cy5.state.calls,file = paste0(path.to.file,'all-particle_Cy5-state-calls.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  utils::write.table(Cy5.dwell.calls,file = paste0(path.to.file,'all-particle_Cy5-dwell-calls.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  utils::write.table(Cy5.residence.calls,file = paste0(path.to.file,'all-particle_Cy5-residence-calls.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  utils::write.table(refined.Cy3.state.calls,file = paste0(path.to.file,'selected_Cy3-state-calls.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  utils::write.table(refined.Cy5.state.calls,file = paste0(path.to.file,'selected_Cy5-state-calls.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  utils::write.table(refined.Cy3.residence.calls,file = paste0(path.to.file,'selected_Cy3-residence-calls.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  utils::write.table(refined.Cy5.residence.calls,file = paste0(path.to.file,'selected_Cy5-residence-calls.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  utils::write.table(refined.Cy3.dwell.calls,file = paste0(path.to.file,'selected_Cy3-dwell-calls.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  utils::write.table(refined.Cy5.dwell.calls,file = paste0(path.to.file,'selected_Cy5-dwell-calls.txt'),quote = FALSE,sep = '\t',col.names = TRUE,row.names = FALSE)
  return(list('residence.data'=residence.data,'Cy3.image.avg'=Cy3.image.avg,'refined.Cy3.traces'=refined.Cy3.traces,'refined.Cy3.trace.rolls'=refined.Cy3.trace.rolls,'REFspots.Cy3axes'=REFspots.Cy3axes,'refined.Cy3.snaps'=refined.Cy3.snaps,'Cy3.state.calls'=Cy3.state.calls,'Cy3.residence.calls'=Cy3.residence.calls,'Cy3.dwell.calls'=Cy3.dwell.calls,'Cy5.image.avg'=Cy5.image.avg,'refined.Cy5.traces'=refined.Cy5.traces,'refined.Cy5.trace.rolls'=refined.Cy5.trace.rolls,'REFspots.Cy5axes'=REFspots.Cy5axes,'refined.Cy5.snaps'=refined.Cy5.snaps,'Cy5.state.calls'=Cy5.state.calls,'Cy5.residence.calls'=Cy5.residence.calls,'Cy5.dwell.calls'=Cy5.dwell.calls))
}





