/*!
\page styleguide Coding Style Guide

Code developed in <B>%bitpit</B> should follow the coding styles described here.  Any deviations from this style
guide will result in severe berating and other verbal abuse.

\section directory_structure bitpit Directory Structure
<B>%bitpit</B> source code is organized in the following directory structure: \n
 - doc: documentation is put here, along with the input file for Doxygen.  Most <B>%bitpit</B> documentation is doxygen-processed;
 - examples: examples of <B>%bitpit</B> usage, both small and large.  These programs are not meant to be used as unit tests, but
     rather as further documentation on <B>%bitpit</B> usage;
 - src : source code;
 - test: all unit test programs should go below this directory.

If you're designing a new class or other code for <B>%bitpit</B> and are not sure where to put it, try to find something similar
and put it there.  Otherwise, email the <B>%bitpit</B> email list for pointers.  

\section source_style Source Code Style and Best Practices

\subsection indentation Indentation
Use for indentation 4 characters.

The preferred way to ease multiple indentation levels in a switch statement is
to align the "switch" and its subordinate "case" labels in the same column
instead of "double-indenting" the "case" labels.  E.g.:
\code{.unparsed}
	switch (suffix) {
	case 'G':
	case 'g':
		mem <<= 30;
		break;
	case 'M':
	case 'm':
		mem <<= 20;
		break;
	case 'K':
	case 'k':
		mem <<= 10;
		fall through
	default:
		break;
	}
\endcode
Don't put multiple statements on a single line unless you have
something to hide:
\code{.unparsed}
	if (condition) do_this;
	  do_something_everytime;
\endcode
Get a decent editor and don't leave whitespace at the end of lines.

\subsection bracing Placing Braces

The other issue that always comes up in C++ styling is the placement of
braces.  Unlike the indent size, there are few technical reasons to
choose one placement strategy over the other, but the preferred way, as
shown to us by the prophets Kernighan and Ritchie, is to put the opening
brace last on the line, and put the closing brace first, thusly:
\code{.unparsed}
	if (x is true) {
		we do y
	}
\endcode
This applies to all non-function statement blocks (if, switch, for,
while, do).  E.g.:
\code{.unparsed}
	switch (action) {
	case KOBJ_ADD:
		return "add";
	case KOBJ_REMOVE:
		return "remove";
	case KOBJ_CHANGE:
		return "change";
	default:
		return NULL;
	}
\endcode
However, there is one special case, namely functions: they have the
opening brace at the beginning of the next line, thus:
\code{.unparsed}
	int function(int x)
	{
		body of function
	}
\endcode
Heretic people all over the world have claimed that this inconsistency
is ...  well ...  inconsistent, but all right-thinking people know that
(a) K&R are _right_ and (b) K&R are right.  Besides, functions are
special anyway (you can't nest them in C++).

Note that the closing brace is empty on a line of its own, _except_ in
the cases where it is followed by a continuation of the same statement,
i.e., a "while" in a do-statement or an "else" in an if-statement, like
this:
\code{.unparsed}
	do {
		body of do-loop
	} while (condition);
\endcode
and
\code{.unparsed}
	if (x == y) {
		..
	} else if (x > y) {
		...
	} else {
		....
	}
\endcode
Rationale: K&R.

Also, note that this brace-placement also minimizes the number of empty
(or almost empty) lines, without any loss of readability.  Thus, as the
supply of new-lines on your screen is not a renewable resource (think
25-line terminal screens here), you have more empty lines to put
comments on.

Do not unnecessarily use braces where a single statement will do.
\code{.unparsed}
	if (condition)
		action();
\endcode
and
\code{.unparsed}
	if (condition)
		do_this();
	else
		do_that();
\endcode
This does not apply if only one branch of a conditional statement is a single
statement; in the latter case use braces in both branches:
\code{.unparsed}
	if (condition) {
		do_this();
		do_that();
	} else {
		otherwise();
	}
\endcode

\subsection naming Naming
Class names should be in the CamelBack style, e.g. PatchCartesian or PatchOctree.

Class methods names should be in CamelBack style style but with a starting lower case letter, e.g. Patch::getId()
or Interface::getArea().

Class private member variables should be in the CamelBack style and should have
the prefix m_, e.g. PatchCartesian::m_cellSize. Each member variable that the
user can modify, should have set/get functions, e.g. for the variable int
m_variable it is necessary to define void setVariable(int newval)
and int getVariable() const.

Enumeration values should be all captitalized, with underscores avoided if possible.
The name of the enumeration name should be in the CamelBack style and indicates
the general purpose of the enumeration, so e.g. we use Interface::PositionType, not
Interface::InterfacePositionType.

No names should be added to the global namespace.  Everything should be
in the <B>%bitpit</B> namespace.  An exception can be made for names with a static
scope declared in a .cpp file, but class member functions never have a
static scope.

Names should be kept as private as possible.  If declaring a struct or
utility class that is used internally by some other class, consider
defining it in the .cpp file of the main class or a separate header
only included in that .cpp file and using (if necessary) only forward
delcarations (e.g. \c class \c Node;) in the header file used
by other code.  If that is not possible, then consider nesting the
definitions such that the scope of the name is limited to that of the
class using it.

Any names introduced into the top-level <B>%bitpit</B> namespace should be
sufficiently unique to avoid conflicts with other code.  If you must
introduce a class to the top-level <B>%bitpit</B> namespace, don't choose
an overly genereric name like \c Node.

\subsection includes Includes

Developers should avoid using \#include in header files, as they propagate
dependencies more widely than necessary.  The only cases where other includes
are needed are to import the declaration for a parent class, and to declare
types used as non-pointer and non-reference function arguments.  In most cases,
a forward-declaration statement (e.g. 'class Interface') will suffice.

\subsection constants-macros Constants and Macros

Don't use a pre-processor macro where a const variable or an inline or
template function will suffice. There is absolutely benefit to the former
over the later with modern compilers.  Further, using  macros bypasses
typechecking that the compiler would otherwise do for you and if used in
headers, introduce names into the global rather than <B>%bitpit</B> namespace.

Don't define constants that are already provided by standard libraries.
For example, use \c M_PI as defined in \c math.h rather than defining
your own constant.

\subsection typedefs Typedefs

Please don't use things like "vps_t".
It's a _mistake_ to use typedef for structures and pointers. When you see a
\code{.unparsed}
	vps_t a;
\endcode
in the source, what does it mean?
In contrast, if it says
\code{.unparsed}
	struct virtual_container *a;
\endcode
you can actually tell what "a" is.

Lots of people think that typedefs "help readability". Not so. They are
useful only for:

 - totally opaque objects (where the typedef is actively used to _hide_
     what the object is).

     Example: opaque objects that you can only access using
     the proper accessor functions.

     <B>NOTE</B> Opaqueness and "accessor functions" are not good in themselves.

 - Clear integer types, where the abstraction _helps_ avoid confusion
     whether it is "int" or "long".

     <B>NOTE</B> Again - there needs to be a _reason_ for this. If something is
     "unsigned long", then there's no reason to do
\code{.unparsed}
	typedef unsigned long myflags_t;
\endcode
     but if there is a clear reason for why it under certain circumstances
     might be an "unsigned int" and under other configurations might be
     "unsigned long", then by all means go ahead and use a typedef.

 - when you use sparse to literally create a _new_ type for
     type-checking.

 - New types which are identical to standard types, in certain
     exceptional circumstances.

     Although it would only take a short amount of time for the eyes and
     brain to become accustomed to the standard types like 'uint32_t',
     some people object to their use anyway.

     When editing existing code which already uses one or the other set
     of types, you should conform to the existing choices in that code.

Maybe there are other cases too, but the rule should basically be to NEVER
EVER use a typedef unless you can clearly match one of those rules.

In general, a pointer, or a struct that has elements that can reasonably
be directly accessed should _never_ be a typedef.

\subsection functions Functions

Functions should be short and sweet, and do just one thing.  They should
fit on one or two screenfuls of text (the ISO/ANSI screen size is 80x24,
as we all know), and do one thing and do that well.

The maximum length of a function is inversely proportional to the
complexity and indentation level of that function.  So, if you have a
conceptually simple function that is just one long (but simple)
case-statement, where you have to do lots of small things for a lot of
different cases, it's OK to have a longer function.

However, if you have a complex function, and you suspect that a
less-than-gifted first-year high-school student might not even
understand what the function is all about, you should adhere to the
maximum limits all the more closely.  Use helper functions with
descriptive names (you can ask the compiler to in-line them if you think
it's performance-critical, and it will probably do a better job of it
than you would have done).

Another measure of the function is the number of local variables.  They
shouldn't exceed 5-10, or you're doing something wrong.  Re-think the
function, and split it into smaller pieces.  A human brain can
generally easily keep track of about 7 different things, anything more
and it gets confused.  You know you're brilliant, but maybe you'd like
to understand what you did 2 weeks from now.

In source files, separate functions with one blank line.

In function prototypes, include parameter names with their data types.
Although this is not required by the C++ language, it is preferred
because it is a simple way to add valuable information for the reader.

\subsection exit-functions Centralized exiting of functions

Albeit deprecated by some people, the equivalent of the goto statement is
used frequently by compilers in form of the unconditional jump instruction.

The goto statement comes in handy when a function exits from multiple
locations and some common work such as cleanup has to be done.  If there is no
cleanup needed then just return directly.

Choose label names which say what the goto does or why the goto exists.  An
example of a good name could be "out_buffer:" if the goto frees "buffer".  Avoid
using GW-BASIC names like "err1:" and "err2:".  Also don't name them after the
goto location like "err_kmalloc_failed:"

The rationale for using gotos is:

- unconditional statements are easier to understand and follow
- nesting is reduced
- errors by not updating individual exit points when making
    modifications are prevented
- saves the compiler work to optimize redundant code away ;)
\code{.unparsed}
	int fun(int a)
	{
		int result = 0;
		char *buffer;

		buffer = kmalloc(SIZE, GFP_KERNEL);
		if (!buffer)
			return -ENOMEM;

		if (condition1) {
			while (loop1) {
				...
			}
			result = 1;
			goto out_buffer;
		}
		...
	out_buffer:
		kfree(buffer);
		return result;
	}
\endcode
A common type of bug to be aware of it "one err bugs" which look like this:
\code{.unparsed}
	err:
		kfree(foo->bar);
		kfree(foo);
		return ret;
\endcode
The bug in this code is that on some exit paths "foo" is NULL.  Normally the
fix for this is to split it up into two error labels "err_bar:" and "err_foo:".


\subsection function-return-values Function return values

Functions can return values of many different kinds, and one of the
most common is a value indicating whether the function succeeded or
failed.  Such a value can be represented as an error-code integer
(-Exxx = failure, 0 = success) or a "succeeded" boolean (0 = failure,
non-zero = success).

Mixing up these two sorts of representations is a fertile source of
difficult-to-find bugs.  To help prevent such bugs, always follow this
convention:

	If the name of a function is an action or an imperative command,
	the function should return an error-code integer.  If the name
	is a predicate, the function should return a "succeeded" boolean.

For example, "add work" is a command, and the add_work() function returns 0
for success or -EBUSY for failure.  In the same way, "PCI device present" is
a predicate, and the pci_dev_present() function returns 1 if it succeeds in
finding a matching device or 0 if it doesn't.

All public functions must respect this convention. Private (static) functions
need not, but it is recommended that they do.

Functions whose return value is the actual result of a computation, rather
than an indication of whether the computation succeeded, are not subject to
this rule.  Generally they indicate failure by returning some out-of-range
result.  Typical examples would be functions that return pointers; they use
NULL to report failure.


\subsection inline The inline disease

There appears to be a common misperception that gcc has a magic "make me
faster" speedup option called "inline". While the use of inlines can be
appropriate (for example as a means of replacing macros), it
very often is not. Abundant use of the inline keyword leads to a much bigger
executable, which in turn slows the program as a whole down, due to a bigger
icache footprint for the CPU and simply because there is less memory
available for the pagecache. Just think about it; a pagecache miss causes a
disk seek, which easily takes 5 milliseconds. There are a LOT of cpu cycles
that can go into these 5 milliseconds.

A reasonable rule of thumb is to not put inline at functions that have more
than 3 lines of code in them. An exception to this rule are the cases where
a parameter is known to be a compiletime constant, and as a result of this
constantness you *know* the compiler will be able to optimize most of your
function away at compile time.

Often people argue that adding inline to functions that are static and used
only once is always a win since there is no space tradeoff. While this is
technically correct, gcc is capable of inlining these automatically without
help, and the maintenance issue of removing the inline when a second user
appears outweighs the potential value of the hint that tells gcc to do
something it would have done anyway.

\subsection conditional-compilation Conditional Compilation

Wherever possible, don't use preprocessor conditionals (\#if, \#ifdef) in .c
files; doing so makes code harder to read and logic harder to follow.  Instead,
use such conditionals in a header file defining functions for use in those .c
files, providing no-op stub versions in the \#else case, and then call those
functions unconditionally from .c files.  The compiler will avoid generating
any code for the stub calls, producing identical results, but the logic will
remain easy to follow.

Prefer to compile out entire functions, rather than portions of functions or
portions of expressions.  Rather than putting an ifdef in an expression, factor
out part or all of the expression into a separate helper function and apply the
conditional to that function.

If you have a function or variable which may potentially go unused in a
particular configuration, and the compiler would warn about its definition
going unused, mark the definition as __maybe_unused rather than wrapping it in
a preprocessor conditional.  (However, if a function or variable *always* goes
unused, delete it.)

At the end of any non-trivial \#if or \#ifdef block (more than a few lines),
place a comment after the \#endif on the same line, noting the conditional
expression used.  For instance:
\code{.unparsed}
	#ifdef CONFIG_SOMETHING
	...
	#endif // CONFIG_SOMETHING
\endcode

\subsection commenting Commenting

Comments are good, but there is also a danger of over-commenting.  NEVER
try to explain HOW your code works in a comment: it's much better to
write the code so that the _working_ is obvious, and it's a waste of
time to explain badly written code.

Generally, you want your comments to tell WHAT your code does, not HOW.
Also, try to avoid putting comments inside a function body: if the
function is so complex that you need to separately comment parts of it,
you should probably split the function into smaller functions.  You can make
small comments to note or warn about something particularly clever (or
ugly), but try to avoid excess.  Instead, put the comments at the head
of the function, telling people what it does, and possibly WHY it does
it.

Try to keep header files free of comments; when comments are inside headers,
all users of those headers must be recompiled if a comment is changed.

Each class should be fully commented with a doxygen comment block. A doxygen
comment block is a special comment block with some additional markings, so
doxygen knows it is a piece of structured text that needs to end up in the
generated documentation. The preferred way to mark a doxygen comment block
is the JavaDoc style, which consist of a C-style comment block starting with
two *'s, like this:
\verbatim
    \**
    * ... text ...
    */
\endverbatim
The doxygen comment block of a class should include a general description of
the class and a list of its features and possible limitations.

A doxygen comment block should be added for every method of a class (both
public and private methos should be commented), this comment block should
contain a general description of the method and a description of all the
arguments and the return value of the method. For instance:
\verbatim
    /**
    * \brief Brief description.
    *        Brief description continued.
    *
    *  Detailed description starts here.
    *
    *  \param var detailed description of the parameter
    *  \result Detailed description of the return value.
    */
\endverbatim

To document the members of a file, struct, union, class, or enum, it is sometimes
desired to place the documentation block after the member instead of before.
For this purpose an additional < marker is required in the comment block.
Note that this also works for the parameters of a function. For instance:
\verbatim
    int var; /**< Detailed description after the member */
\endverbatim
To simplify the documentation of class members, define just one member per
line (no commas for multiple members declarations).

Doxygen commands should start with a backslash (\\).

As a rule of thumb, the code should run through doxygen without generating any
warnings; in fact, doxygen is sometimes helpful at pointing out inconsistencies
in your class declaration.

\section git Git Repository Practices
As most of our code repositories uses git as the revision control system, it is important to decide on a workflow that can be followed by the individual developer. The way that any individual developer interact with the upstream git repository can have an important impact on other developers and the ability to identify and manage individual changes.  This set of guidelines and practices attempts to establish some standards for how developers will interact with the upstream git repository.

\subsection commits Making Repository Commits
As a general rule, developers should update frequently, and commit changes often.  However, the repository should always remain
in a state where the code can be compiled.  Most of the time, the code should also successfully execute "make check" run from the
top-level directory.  If you commit code that violates this principal, it should be your first priority to return the repository
code to a compilable state, and your second priority to make sure "make check" runs without errors.
Although it would be possible and many software projects do it, we prefer not to force successful execution of the test suite
before every commit.  Developers should make every effort to avoid having to impose this constraint, by running a make check
before every commit.

Commits to the repository should also come with a non-trivial and useful log message.

\subsection outside-master Working Outside the Master Branch
A critical concept is that all changes shall be developed outside of the master<sup>1</sup> branch.  Whether they are in a different branch of the upstream<sup>2</sup> repository (gitflow) or a branch of an entirely different fork (forking workflow) is secondary.  This is a well-established concept regardless of the workflow being adopted, and allows a number of other benefits as described below.

\subsubsection fork Working on a Different Fork
There are a number of benefits of working on a different fork rather than a branch of the upstream repo, although not strictly technical:
- Developers, particularly new developers, are liberated from the constant oversight of others as they explore new code options.  The impact of this may depend on an individual developer’s personality, but for some it creates a refuge where they can be more free and creative.
- Similarly, assuming that all changesets in the upstream repo are communicated to the entire development team, the team is spared a noisy stream of notifications and can focus their attention on the rarer occurrence of a pull request notification.

\subsubsection pr All Changes are Committed by Pull Request
Although this can be imposed technically by limiting the authority to alter the upstream repo (as in the forking workflow), a healthy developer community can also simply rely on convention.  The advantage of doing it by convention rather than by restriction is that it is easier to distribute the load of reviewing and accepting changes.
A critical consequence of this decision is that all code is reviewed before it is committed to the upstream master branch.  This has benefits to overall quality in two related ways:
- the code under review will improve due to the review itself, and
- those involved in the review will maintain a broad awareness of the code base resulting in better contributions from them.

This practice does, however, place a substantial burden on the developers to perform timely reviews of the pull requested (PR’ed) code.  PR’s that languish without sufficient review have a number of negative consequences:
- they need to be refreshed simply to keep them up-to-date with the possibly advancing upstream/master
- they may delay further development on similar or related features
- they breed frustration in the original developer, undermining the community as a whole.
github provides powerful collaboration tools that greatly facilitate this process.

<sup>1</sup> Although a repository may choose a different name for its main development branch, this document will refer to that as the “master” branch.

<sup>2</sup> For this discussion, the “upstream” repo will refer to the centralized authoritative repository used to synchronize changes.

\subsection git-mechanics Some Git Mechanics to Keep it Clean
Given the above practices, there are some mechanical details that can help ensure that the upstream/master repository is always in a state that facilitates all repository actions and interactions.

-# Feature branches being used for development should be kept up-to-date with the upstream/master by rebase only.  When a feature branch is rebased against the upstream/master, all changes in the upstream/master are inserted into the feature branch at a point in its history that is prior to any of the changes of the feature branch.  This can require conflict resultion as the feature branch changes are “replayed” on top of the new upstream/master in its current state.  The primary advantage of this policy is that it keeps all of the feature branch changes contiguous.  If, by contrast, the upstream/master is merged into the feature branch, the recent changes in the upstream/master become woven into the set of changes in the feature branch.  This can make it more difficult to isolate changes later on.

    Strict adoption of this practice is important since a single merge into a feature branch that is then merged back into the upstream/master can make it nearly impossible for others to rebase.

    A typical workflow with pull-request might look like this, all using the command-line, except for submitting the final pull request.  Note that there is never a merge operation.
    -# synchronize your local `master` branch before anything else
    \code{.sh}
     $ git checkout master
     $ git fetch upstream
     $ git rebase upstream/master
     \endcode

    -# now create a new feature branch from master
    \code{.sh}
     $ git checkout -b my_feature_branch master
    \endcode

    -# now make changes, editing A.cpp, B.hpp, C.cpp

    -# now add/commit your changes to your local feature branch
    \code{.sh}
     $ git add A.cpp B.hpp C.cpp
     $ git commit -m “Make sure you have a good commit message”
    \endcode
    -# push your changes to your feature branch on your fork (often called `origin`)
    \code{.sh}
     $ git push origin my_feature_branch
    \endcode
    -# make more changes, editing B.hpp, D.hpp, E.cpp

    -# add/commit your changes to your local feature branch
    \code{.sh}
    $ git add B.hpp D.hpp E.cpp
    $ git commit -m “Be sure you have another good commit message”
    \endcode
    -# push your changes to your freature branch on your fork (often called `origin`)
    \code{.sh}
    $ git push origin my_feature_ranch
    \endcode
    -# When you are ready to submit a pull request, be sure that your feature branch is up-to-date. This first step may seem redundant but is here to be clear which branch we are acting on
    \code{.sh}
    $ git checkout my_feature_branch
    $ git fetch upstream
    $ git rebase upstream/master
    \endcode
      This may generate conflicts that can be addressed at this point.

      NOTE: This step can be performed at any time and should be performed as often as practical to reduce the scope of potential conflicts.

    -# push your updated feature branch on your fork (often called `origin`)
     \code{.sh}
     $ git push origin my_feature_branch
     \endcode
      This may require the ‘-f’ option to force the push.  (It is frequently necessary to force this push because the act of rebasing will “replay” the commits from the feature branch on top of the master, leading to different commit hashes.  Each of the commits will contain the same actual information, but because it has a different set of hashes, git will think there is an inconsistency and ask you to force the change.)

    -# Submit a pull request on github, from your fork to the fathomteam fork.


-# When ready to be adopted into the upstream/master, feature branches should be combined by merge only.  This adds the changeset to the end of the upstream/master as a set of individual commits but in a contiguous block.

   A typical workflow to merge a pull-request might look like this, all using the command-line.
   -# synchronize your local `master` branch before anything else (just because it’s never a bad idea!)
   \code{.sh}
   $ git checkout master
   $ git fetch upstream
   $ git rebase upstream/master
   \endcode
   -# add a remote for the user with the pull-request, perhaps the user is ‘other_user’
   \code{.sh}
   $ git remote add other_user \
         git@bitbucket.org:other_user/moab.git
   \endcode
   -# fetch the other users repo
   \code{.sh}
   $ git fetch other_user
   \endcode
   -# check out their feature branch
   \code{.sh}
   $ git checkout -b pr_feature_branch \
       other_user/feature_branch
   \endcode
   -# confirm that it is up-to-date with the master. This first step may seem redundant but is here to be clear which branch we are acting on
   \code{.sh}
   $ git checkout pr_feature_branch
   $ git fetch upstream
   $ git rebase upstream/master
   \endcode
   This may generate conflicts that can be addressed at this point.  You may want to request the original author (other_user) take care of these.
   -# once confirmed that it’s up-to-date with master, review this branch including:
      -reading the code
      -building the code
      -running tests
   -# once satisfied that the code meets the necessary standards and that all required/requested changes are fully incorporated into other_users’s feature branch, merge it into master
    \code{.sh}
    $ git checkout master
    \endcode
    The next two steps may seem redundant but provide some QA
    \code{.sh}
    $ git fetch upstream
    $ git rebase upstream/master
    $ git merge other_user/feature_branch
    \endcode
   -# push those changes into the master branch on bitbucket
    \code{.sh}
    $ git push upstream/master
    \endcode

-# When a pull request is open for review, any changes to the feature branch will automatically update the pull request.  This is the appropriate way for a developer to respond to requests for changes that occur through the PR process.


-# If a developer has ongoing work that is based on a feature branch that is under consideration in an open PR, a new feature branch (B) should be created that is based on the previous feature branch (A).  Moreover, as changes are made to the original feature branch (A) due to the review process, the new feature branch (B) should be kept up-to-date by rebase against feature branch (A).  This keeps all subsequent changes of (B) downstream from the changes in (A).  Once feature branch (A) has been adopted into the upstream/master, the new feature branch (B) can start being rebased against the upstream/master instead.


-# When a repo is forked, its branches are not automatically synchronized with the corresponding branches on the upstream repo.  This requires a manual process of synchronization via a local clone.  Assuming that the local repo’s branch has the same name as the upstream branch (\<branch\>), and that the fork is known as “origin”:
    \code{.sh}
     $ git fetch upstream
     $ git checkout <branch>
     $ git rebase upstream/<branch>
     $ git push origin <branch>
    \endcode
The decision of which branches to keep up-to-date is up to the developers.  Developers may choose to delete some branches from their own fork to avoid (a) the need to update it and (b) accidentally assuming that it is up-to-date.


-# When rebasing, it is not uncommon to encounter conflicts.  This will interrupt the rebasing process, and each conflicted file will contain conflict blocks.  You will be offered three choices:
 - manually resolve the conflicts, add the resolved files with git add,  and then git rebase --continue (do not commit the resolved files!)
 - abort the rebasing process entirely with git rebase --abort
 - skip the commit that causes the conflict, assuming that you are sure that this is the right thing to do, with git rebase --skip


-# github offers a number of buttons/tools to manage changes between branches and forks.  Some of these operate in ways that are contrary to the practices recommended here, and others are consistent with these practices.  In general, it is best to know how to do all of those operations with the command-line instead of relying on github, as it gives you full control over important details.


-# During code development, it might be necessary to work on the same branch on different machines. The workflow to update the local branch is to first fetch the remote changes and then perform a hard reset.
     \code{.sh}
     $ git fetch origin
     $ git reset --hard origin/branch_name
     \endcode
One should be careful with the branch name as a hard reset would overwrite all changes in the working directory.


Top: \ref styleguide

*/
