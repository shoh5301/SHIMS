#include "Calc.h"

double stack[MAX];
int top;

int inner_product(int v1[],int v2[]){
	return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

double eval(char* post){
//    printf("%s", post);
    char t[MAX];   
    precedence token;
    double temp;
    int i=0,j=0,cnt_operand=0,cnt_operator=0;
    init_stack();
    while((token = get_token(post))!=eos){
        if(token == space && i>0){ //공백을 만나면 token을 실수형으로 변환
            t[i] =  '\0';
            temp = atof(t);
            push(temp);
            cnt_operand++;
            memset(t,0,MAX);
            i=0;
        }//문자가 아니고 operand이거나 '-'이면서 다음에 operand이면 실행
        else if((!isalpha(*post) && operand == token) || (token == minus && get_token(post+1) == operand))
            t[i++] = *post;
        else{   //operator이면 switch문 실행
            if(!(token >= Sin && token <= space) && !isalpha(*post) && token != lparen)
                cnt_operator++;
            switch(token){
                case Sin    : push(sin(pop()*M_PI/180));
                              break;
                case Cos    : push(cos(pop()*M_PI/180));
                              break;
                case Tan    : push(tan(pop()*M_PI/180));
                              break;
                case ln     : push(log(pop()));
                              break;
                case plus   : push(pop() + pop());
                              break;
                case minus  : temp = pop();
                              push(pop() - temp);
                              break;
                case square:  temp = pop();
                              push(pow(pop(), temp));
                              break;
                case times  : push(pop() * pop());
                              break;
                case divide : if((temp = pop()) != 0);
                              push(pop() / temp);
                              break;
                case mod    : temp = pop();
                              push((int)pop() % (int)temp);
                              break;
            }
        }
        post++;
    }

    if(cnt_operand - cnt_operator != 1){
        puts("ERROR: Wrong parameter format...");
        exit(1);
    }
    return pop();
}

void change_postfix(char* str, char* postfix){
    int i=0,j=0,flag=0,count=0,s_c_t;
    precedence token;
  int icp[] = {20, 19, 12, 12, 13, 13, 13, 18, 0, 20, 20, 20, 20}; //스택에 들어가는 값의 우선순위
  int isp[] = {0, 19, 12, 12, 13, 13, 13, 18, 0, 0, 0, 0, 0};     //스택에 있는 값의 우선순위

	init_stack();
    push(eos);    //우선순위가 가장 낮은 '\0'를 삽입
    while((token = get_token(str + i)) != eos){ // str이 널값이면 변환 중지
        if(token == operand){ // operand이면 postfix에 바로 삽입
            if(s_c_t = sin_cos_tan(str + i)){
                if(count%2){  //-1 = (0-1) 이기때문에 ==>  0 1 - 로 표현
                    strncat(postfix+j, "0 ", 2);  //operand 뒤에 공백 삽입
                    j +=2;
                    push(minus);
                    //push(token);
                    count = 0;
                }
                strncat(postfix+j, s_c_t == Sin ? "sin" : s_c_t == Cos ? "cos" : "tan", 3);
                push(s_c_t);
                if(lparen == get_token(str + i + 3)){
                    strncat(postfix+j, "( ", 2);
                    i++;
                    j +=2;
                }
                j += 3;
                i += 3;
                continue;
            }
            else if (is_ln(str+i)) {
                if (count % 2) {  //-1 = (0-1) 이기때문에 ==>  0 1 - 로 표현
                    strncat(postfix + j, "0 ", 2);  //operand 뒤에 공백 삽입
                    j += 2;
                    push(minus);
                    //push(token);
                    count = 0;
                }
                strncat(postfix + j, "ln", 2);
                push(ln);
                if (lparen == get_token(str + i + 2)) {
                    strncat(postfix + j, "( ", 2);
                    i++;
                    j += 2;
                }
                j += 2;
                i += 2;
                continue;
            
            }
            else if (count && count % 2) {
                *(postfix + j++) = '-';
                count = 0;
            }
            *(postfix + j++) = *(str + i);
            flag = 1;     //operand가 실행 되었다는것을 표시
        }else if(token == lparen && count%2){
            strncat(postfix+j, "0 ", 2);  //operand 뒤에 공백 삽입
            j +=2;
            push(minus);
            push(token);
            count = 0;
        }
        //음수 판별 조건 : token은 '-' 이고, operand 또는(or)     ')' 바로 뒤에 '-'는 뺄셈 부호임
        //                 따라서 operand가 선행 안되면서(and),   ')'가 이전에 없으면 '-'는 부호 표시 기호이다.
        else if(token == minus && !flag && get_token(str + i - 1) != rparen)
            count++;
        else{
            if(flag){ //operand와 구분 시켜 주기위해 공백 삽입, flag가 1이면 operand가 삽입 되었단 말임
                *(postfix + j++) = ' ';
                flag = 0;
            }
            if(token == rparen){ // token이 ')'이면 '('이 나올때까지 postfix에 삽입
                while((token = (int)pop()) != lparen && !(token >= Sin && token <= ln)){
                    *(postfix + j++) = ' ';
                    *(postfix + j++) = print_token(token);
                    *(postfix + j++) = ' ';
                }
                if(token != lparen){  //sin => $ , cos => @ , tan => # 기호 삽입
                    *(postfix + j++) = print_token(token);
                    *(postfix + j++) = ' ';
                }
            }else{ // token값이 stack[top]값보다 우선순위가 높을때까지 postfix에 삽입
                while(icp[token] <= isp[(int)stack[top]]){
                    *(postfix + j++) = ' ';
                    *(postfix + j++) = print_token((int)pop());
                    *(postfix + j++) = ' ';
                }
                push(token); //마지막에 token값을 푸쉬
            }
        }
        i++;
    }
   
    *(postfix + j++) = ' ';  //operand와 구분 시켜 주기위해 공백 삽입 
   
    while ((token = (int)pop()) != eos)   //스택에 남은 값들을 postfix에 삽입
    {
        *(postfix + j++) = print_token(token);
    }
    *(postfix + j) = '\0'; //postfix 맨마지막에 null 삽입

    return;
}

precedence is_ln(char* ptr) {
    if (!strncmp(ptr, "ln", 2))
        return ln;
    return 0;
}

precedence sin_cos_tan(char* ptr){
    if(!strncmp(ptr,"sin",3))
        return Sin;
    if(!strncmp(ptr,"cos",3))
        return Cos;
    if(!strncmp(ptr,"tan",3))
        return Tan;
    return 0;
}

precedence get_token(char* symbol){ //token
    switch(*symbol){
        case '(' : return lparen;
        case ')' : return rparen;
        case '+' : return plus;
        case '-' : return minus;
        case '*' : return times;
        case '/' : return divide;
        case '%' : return mod;
        case '\0': return eos;
        case ' ' : return space;
        case '$' : return Sin;
        case '@' : return Cos;
        case '#' : return Tan;
        case '^' : return square;
        case '!' : return ln;
        default  : return operand;
    }
}

char print_token(precedence token){
    switch(token){
        case plus : return '+';
        case minus : return '-';
        case times : return '*';
        case divide : return '/';
        case mod : return '%';
        case Sin : return '$';
        case Cos : return '@';
        case Tan : return '#';
        case ln: return '!';
        case square: return '^';
        default  : return '\0';
    }
}

void init_stack(void){
    top = -1;
}

void push(double in){
    if(top >= MAX-1){
        puts("Stack is Full !!!!!!!!!!!");
        exit(1);
    }
    stack[++top] = in;
}

double pop(void){
    if(top < 0){
        puts("Stack is Empty !!!!!!!!!!!");
        exit(1);
    }
    return stack[top--];
}

void Etodec(char* eq){ //Change 3E4 -> 3*10^4
	int i,prev,post;
	char* ptr;
	char temp[MAX]={'\0'},num[10];
	double dec;

	sprintf(num,"*10^");

  while(1){
  	ptr=NULL;
  	ptr=strstr(eq,"E");
	if(ptr==NULL)
		return;
	else{
		for(i=0;i<MAX;i++){
			if(eq[i]=='E')
				break;
		}
		prev=i;
		strcpy(temp,eq);
		for(i=prev;num[i-prev]!='\0';i++)
			eq[i]=num[i-prev];
		if(temp[i-3]=='+'){
		    for(i=prev+4;temp[i-2]!='\0';i++) //get chars after exponent
			eq[i]=temp[i-2];
		    eq[i]='\0';
		}else{
		    for(i=prev+4;temp[i-3]!='\0';i++) //get chars after exponent
			eq[i]=temp[i-3];
		    eq[i]='\0';
		}
	}
  }
	return;
}

void time_ins(char* Gline,double const time){
        char temp[MAX]={'\0'},temp2[MAX]={'\0'};
        char* tptr;
        int i;

        while(1){
                tptr=strstr(Gline,"t");
                if(tptr!=NULL){
                        strcpy(temp,tptr);
                        for(i=1;i<MAX;i++){
                                temp[i-1]=temp[i];
                                if(temp[i]=='\0')
                                        break;
                        }
                        sprintf(temp2,"%.2f",time);
                        strcpy(tptr,temp2);
                        strcat(Gline,temp);
                }else
                        break;
        }
        return;
}       

