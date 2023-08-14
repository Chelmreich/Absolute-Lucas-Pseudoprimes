class Stack
{
  int items[20];
  int top;

  public:
  Stack() { top=0; }
  inline void push(int x) { items[top++]=x; }
  inline int pop() { return items[--top]; }
  inline bool isempty() { return(top==0); }
  inline int gettop() { return items[top-1]; }
};

