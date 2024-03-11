# Introduction
The purpose of this style guide is to describe stylistic conventions and to ensure consistency and uniformity throughout the code. In general, we will follow the rules established in Google's Python Style Guide. Additionally, there are the SUGAR-specific rules that should be followed, and that supersede Google Python Style Guide. The remaining sections highlight the SUGAR-specific guidelines.
 
## Variables
 Nodes and state variables should adhere to the practices outlined below. All other variable names should be as descriptive as possible. Strive for clarity over brevity. If you must use an acronym or shorthand because you find that a variable name is too long (25 characters or more), then feel free to do so. If the shorthand or acronym is not apparent, define it in a comment.

### State Variables
Below are the variable names that should only be used to define state variables.
* V - Voltage
* I - Current
* P - Active Power
* Q - Reactive Power
* G - Conductance
* B - Susceptance
* X - Reactance
* R - Resistance
* Y - Admittance
* Z - Impedance

If a state variable is phase-specific, define it like this:
```
Ib = complex(10, 11) # phase B current expressed in complex plane
In = complex(0, 0)  # neutral phase current expressed in complex plane
V1 = complex(120, 0) # phase 1 triplex voltage
```

If a state variable has a real or imaginary component and a phase component, define it like this:
```
Var = 1.0  # real phase A voltage
Vbr_from = 2.0 # the phase B sending voltage
Vbr_to = 2.0 # the phase B receiving voltage
```

If a state variable is defined in terms of polar coordinates, use the format below:
```
Ib_mag = 5.0 # magnitude of the phase B current
Ib_ang = 10.0 # angle in degrees of the phase B current
```

### Nodes
 Node names should be formatted as follows:
    `node{phase}_{state variable}{r/i}{p/s}_{to/from or descriptor}`
 Here, `{Phase}` refers to whether the node is a three-phase (A, B or C) or triplex (1, 2 or N) node. `r/i` indicates if
  the node
  corresponds to
 the real or imaginary component of a variable.`p/s` is used when the variable is on a primary or secondary winding. `to
 /from` is used if the model the node is within is
 bidirectional. If the node is model specific, use a `descriptor` that clarifies the purpose of the node, such as the
  element name or the node's purpose. If the state variable is inherently real or
  imaginary (e.g., P and Q), then there is no need to adhere to the `r/i` convention. 

 Nodes that apply to all models are created in [Nodes.py](../classes/Nodes.py). Nodes that are class-specific should
  be created within their respective class within a `def assign_nodes()` function.

Here are examples of excellent and lousy node names:
```python

Good:        
    nodeA_Vr_from
    nodeB_Ir_to
    nodeA_Vr_fuse
    nodeC_Ir_measure
    nodeGnd_Vrp
    node1_Vr_from
Bad:
    node_fuseA_Vr
    nodeXtraA_rs
    nodeVGnd_rp
    nodeVr_A_from     
```

## Functions
All functions and methods should be lowercase and use underscores for clarity. Also, functions should start with a
 docstring that provides a summary, a description of function arguments, and specify any function returns. The examples provided by the Google Python Style Guide are below. Either format is acceptable.
  
  ```python
def fetch_smalltable_rows(table_handle: smalltable.Table,
                          keys: Sequence[Union[bytes, str]],
                          require_all_keys: bool = False,
                         ) -> Mapping[bytes, Tuple[str]]:
    """Fetches rows from a Smalltable.

    Retrieves rows pertaining to the given keys from the Table instance represented by table_handle.  String keys will be UTF-8 encoded.

    Args:
        table_handle: An open smalltable.Table instance.
        keys: A sequence of strings representing the key of each table row to fetch.  String keys will be UTF-8 encoded.
        require_all_keys: Optional; If require_all_keys is True only rows with values set for all keys will be returned.

    Returns:
        A dict mapping keys to the corresponding table row data
        fetched. Each row is represented as a tuple of strings. For
        example:

        {b'Serak': ('Rigel VII', 'Preparer'),
         b'Zim': ('Irk', 'Invader'),
         b'Lrrr': ('Omicron Persei 8', 'Emperor')}

        Returned keys are always bytes.  If a key from the keys argument is
        missing from the dictionary, then that row was not found in the
        table (and require_all_keys must have been False).

    Raises:
        IOError: An error occurred accessing the smalltable.
    """
```
```python
def fetch_smalltable_rows(table_handle: smalltable.Table,
                          keys: Sequence[Union[bytes, str]],
                          require_all_keys: bool = False,
                         ) -> Mapping[bytes, Tuple[str]]:
    """Fetches rows from a Smalltable.

    Retrieves rows pertaining to the given keys from the Table instance represented by table_handle.  String keys will be UTF-8 encoded.

    Args:
      table_handle:
        An open smalltable.Table instance.
      keys:
        A sequence of strings representing the key of each table row to fetch.  String keys will be UTF-8 encoded.
      require_all_keys:
        Optional; If require_all_keys is True only rows with values set for all keys will be returned.

    Returns:
      A dict mapping keys to the corresponding table row data
      fetched. Each row is represented as a tuple of strings. For
      example:

      {b'Serak': ('Rigel VII', 'Preparer'),
       b'Zim': ('Irk', 'Invader'),
       b'Lrrr': ('Omicron Persei 8', 'Emperor')}

      Returned keys are always bytes.  If a key from the keys argument is
      missing from the dictionary, then that row was not found in the
      table (and require_all_keys must have been False).

    Raises:
      IOError: An error occurred accessing the smalltable.
    """
```
  

### Stamping
 Stamping function names should always be `stamp_linear` for linear power flow elements or for nonlinear elements with
  a linear stamp and `stamp_nonlinear` for nonlinear elements.

## Classes
All class names should use pascal scale and begin with a docstring that provides a summary, a longer description, and
 the attributes of the class.
```python
class NewClass():
    """Summary of class here.

    Longer class information....
    Longer class information....

    Attributes:
        is_cool: A boolean indicating if someone is cool or not.
        cool_factor: An integer count of someone's cool rating.
    """

    def __init__(self, is_cool=False):
        """Inits SampleClass with blah."""
        self.is_cool = is_cool
        self.cool_factor = 0

    def public_method(self):
        """Performs operation blah.""" 
```
Linear models should inherit the [LinearElement](../classes/Elements.py) class.
```python
class NewClass(LinearElement):
    """Summary of class here.

    Longer class information....
    Longer class information....

    Attributes:
        is_cool: A boolean indicating if someone is cool or not.
        cool_factor: An integer count of someone's cool rating.
    """

    def __init__(self, is_cool=False):
        """Inits SampleClass with blah."""
        super(NewClass, self).__init__()
        self.is_cool = is_cool
        self.cool_factor = 0

    def public_method(self):
        """Performs operation blah."""  
```
Nonlinear models should inherit the [NonlinearElement](../classes/Elements.py) class.
```python
class NewClass(NonlinearElement): 
        """Summary of class here.

    Longer class information....
    Longer class information....

    Attributes:
        is_cool: A boolean indicating if someone is cool or not.
        cool_factor: An integer count of someone's cool rating.
    """

    def __init__(self, is_cool=False):
        """Inits SampleClass with blah."""
        super(NewClass, self).__init__()
        self.is_cool = is_cool
        self.cool_factor = 0

    def public_method(self):
        """Performs operation blah."""  
    pass
```

## Comments
If it's not obvious, add a comment. Only use docstrings for file, class, or function descriptions.
### TODO Comments
 Guidelines on TODO comments should adhere to the Google Python Style Guide. The guidelines are reprinted below to underscore their importance.
 
Use `TODO` comments for temporary code, a short-term solution, or
good-enough but not perfect.

A `TODO` comment begins with the string `TODO` in all caps and a parenthesized
name, email address, or another identifier
of the person or issue with the best context about the problem. This is followed
by an explanation of what there is to do.

The purpose is to have a consistent `TODO` format that can be searched to find
out how to get more details. A `TODO` is not a commitment that the person
referenced will fix the problem. Thus when you create a
`TODO`, it is almost always your name
that is given.

```python
# TODO(kl@gmail.com): Use a "*" here for string repetition.
# TODO(Zeke) Change this to use relations.
```

If your `TODO` is of the form "At a future date, do something" make sure that you
either include a precise date ("Fix by November 2009") or a very specific
event ("Remove this code when all clients can handle XML responses.").

## Files
File names should be descriptive but not overly long. A file whose main contents are a `class` or a `function` should have to use the title of the class or function as their name.
 
 All python files (except for runSUGAR3.py) should begin with a docstring. The format for docstrings is:
 ```python
"""A one-line summary of the module or program, terminated by a period.

Author(s): John Doe, Jane Doe
Created Date: 10-10-2019
Updated Date: 06-20-2020
Email: johndoe@cmu.edu, janedoe@cmu.edu
Status: Development  
# other options are: Validated, Deprecated. Use deprecated when something needs an overhaul or
# is outdated, validated if the program has passed validation protocol, and development if it is a work-in-progress.

Leave one blank line.  The rest of this docstring should contain an
overall description of the module or program.  Optionally, it may also
include a brief description of exported classes and functions and/or usage
examples.

  Typical usage example:

  foo = ClassFoo()
  bar = foo.FunctionBar()
"""
```

 If you've contributed anything that is not a minor bug fix or formatting correction, you should add your name to the
 author list.

## Final Words
 You should refer to Google's Python Style Guide for any items not mentioned here such as line length (limit is 80)
 , indentation, and import formatting. Google's YAPF formatter is installed in this repository and can be run at any
  time from the project's root folder using the `make reformat` command. It will ensure that all of the Google Python
   Style Guide rules are enforced.
  
  If you see code in this repository that does not follow these guidelines, feel free to fix it and bring it up to standard.
