function y = tone(x,j)
  switch x
        case 0
       y= sin(2*pi*941*j)+sin(2*pi*1336*j);
            case 1
         y= sin(2*pi*697*j)+sin(2*pi*1209*j);
          case 2
         y= sin(2*pi*697*j)+sin(2*pi*1336*j);    
          case 3
         y= sin(2*pi*697*j)+sin(2*pi*1477*j);    
          case 4
         y= sin(2*pi*770*j)+sin(2*pi*1209*j);    
          case 5
         y= sin(2*pi*770*j)+sin(2*pi*1336*j);    
          case 6
         y= sin(2*pi*770*j)+sin(2*pi*1477*j);    
          case 7
         y= sin(2*pi*852*j)+sin(2*pi*1209*j);    
          case 8
         y= sin(2*pi*852*j)+sin(2*pi*1336*j);    
          case 9
         y= sin(2*pi*852*j)+sin(2*pi*1477*j);    
         
        end     

end